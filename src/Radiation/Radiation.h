// ----------------------------------------------------------------------------
//! \file Radiation.h
//
//! \brief This file contains the header for the generic class Radiation
//   for the particle radiation losses.
//
// ----------------------------------------------------------------------------

#ifndef RADIATION_H
#define RADIATION_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "Params.h"
#include "Particles.h"
#include "Species.h"
#include "RadiationTables.h"
#include "Random.h"

//  ----------------------------------------------------------------------------
//! Class Radiation
//  ----------------------------------------------------------------------------
class Radiation
{

public:
    //! Creator for Radiation
    Radiation( Params &params, Species *species, Random * rand );
    virtual ~Radiation();

    //! Overloading of () operator
    //! \param particles   particle object containing the particle
    //!                    properties of the current species
    //! \param photons     Particles object that will received the emitted photons
    //! \param smpi        MPI properties
    //! \param RadiationTables Cross-section data tables and useful functions
    //                     for nonlinear inverse Compton scattering
    //! \param istart      Index of the first particle
    //! \param iend        Index of the last particle
    //! \param ithread     Thread index
    //! \param radiated_energy     overall energy radiated during the call to this method
    virtual void operator()(
        Particles       &particles,
        Particles       *photons,
        SmileiMPI       *smpi,
        RadiationTables &RadiationTables,
        double          &radiated_energy,
        int             istart,
        int             iend,
        int             ithread,
        int             ibin = 0,
        int             ipart_ref = 0) = 0;

    //! Computation of the Lorentz invariant quantum parameter
    //! for the given particle properties
    //! \param charge_over_mass2 charge divided by the square of the mass
    //! \param px particle x momentum
    //! \param py particle y momentum
    //! \param pz particle z momentum
    //! \param gamma particle Lorentz factor
    //! \param Ex x component of the particle electric field
    //! \param Ey y component of the particle electric field
    //! \param Ez z component of the particle electric field
    //! \param Bx x component of the particle magnetic field
    //! \param By y component of the particle magnetic field
    //! \param Bz z component of the particle magnetic field
    //#pragma omp declare simd
    inline double __attribute__((always_inline)) computeParticleChi( double charge_over_mass2,
                                      double px, double py, double pz,
                                      double gamma,
                                      double Ex, double Ey, double Ez,
                                      double Bx, double By, double Bz )
    {

        return std::fabs( charge_over_mass2 )*inv_norm_E_Schwinger_
               * std::sqrt( std::fabs(  (Ex*px + Ey*py + Ez*pz) * (Ex*px + Ey*py + Ez*pz)
                             - (gamma*Ex - By*pz + Bz*py) * (gamma*Ex - By*pz + Bz*py)
                             - (gamma*Ey - Bz*px + Bx*pz) * (gamma*Ey - Bz*px + Bx*pz)
                             - (gamma*Ez - Bx*py + By*px) * (gamma*Ez - Bx*py + By*px) ) );
    };

    //! Computation of the quantum parameter for the given
    //! thread of particles
    //! \param Particles class containg the particle property arrays
    //! \param smpi class for mpi parameters
    //! \param istart      Index of the first particle
    //! \param iend        Index of the last particle
    //! \param ithread     Thread index
    //#pragma acc routine seq
    void computeParticlesChi( Particles &particles,
                              SmileiMPI *smpi,
                              int istart,
                              int iend,
                              int ithread,
                              int ipart_ref = 0 );
                              
protected:

    // ________________________________________
    // General parameters

    //! Dimension of position
    int n_dimensions_;

    //! Inversed species mass
    double one_over_mass_;

    //! Time step
    double dt_;

    Random * rand_;

    // _________________________________________
    // Factors

    //! Normalized Schwinger Electric field
    double norm_E_Schwinger_;

    //! Inversed Normalized Schwinger Electric field
    double inv_norm_E_Schwinger_;

    //! Particle dimension
    int nDim_;

private:

};//END class

#endif
