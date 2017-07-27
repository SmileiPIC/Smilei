// ----------------------------------------------------------------------------
//! \file MultiphotonBreitWheeler.h
//
//! \brief This file contains the class functions for the generic class
//!  MultiphotonBreitWheeler for the photon decay into pairs via the mutliphoton
//!  Breit-Wheeler process.
//
// ----------------------------------------------------------------------------

#ifndef MULTIPHOTONBREITWHEELER_H
#define MULTIPHOTONBREITWHEELER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "MultiphotonBreitWheelerTables.h"

//  ----------------------------------------------------------------------------
//! Class Radiation
//  ----------------------------------------------------------------------------
class MultiphotonBreitWheeler
{
    public:

        //! Creator for Radiation
        MultiphotonBreitWheeler(Params& params, Species *species);
        ~MultiphotonBreitWheeler();

        //! Overloading of () operator
        //! \param particles   particle object containing the particle
        //!                    properties of the current species
        //! \param smpi        MPI properties
        //! \param MultiphotonBreitWheelerTables Cross-section data tables and useful functions
        //                     for multiphoton Breit-Wheeler
        //! \param istart      Index of the first particle
        //! \param iend        Index of the last particle
        //! \param ithread     Thread index
        void operator() (
                Particles &particles,
                SmileiMPI* smpi,
                MultiphotonBreitWheelerTables &MultiphotonBreitWheelerTables,
                int istart,
                int iend,
                int ithread);

        //! Computation of the photon Lorentz invariant quantum parameter
        //! for the given photon properties
        //! \param kx photon x momentum
        //! \param ky photon y momentum
        //! \param kz photon z momentum
        //! \param gamma photon Lorentz factor
        //! \param Ex x component of the particle electric field
        //! \param Ey y component of the particle electric field
        //! \param Ez z component of the particle electric field
        //! \param Bx x component of the particle magnetic field
        //! \param By y component of the particle magnetic field
        //! \param Bz z component of the particle magnetic field
        //#pragma omp declare simd
        double inline compute_chiph(double & kx, double & ky, double & kz,
                                    double & gamma,
                                    double & Ex, double & Ey, double & Ez,
                                    double & Bx, double & By, double & Bz)
        {

            return inv_norm_E_Schwinger
                  * sqrt( fabs( pow(Ex*kx + Ey*ky + Ez*kz,2)
                  - pow(gamma*Ex - By*kz + Bz*ky,2)
                  - pow(gamma*Ey - Bz*kx + Bx*kz,2)
                  - pow(gamma*Ez - Bx*ky + By*kx,2)));
        };

        //! Computation of the quantum parameter for the given
        //! thread of photons
        //! \param Particles class containg the particle property arrays
        //! \param smpi class for mpi parameters
        //! \param istart      Index of the first particle
        //! \param iend        Index of the last particle
        //! \param ithread     Thread index
        void compute_thread_chiph(Particles &particles,
                SmileiMPI* smpi,
                int istart,
                int iend,
                int ithread);

    private:

        // ________________________________________
        // General parameters

        //! Dimension of position
        int nDim_;

        //! Time step
        double dt;

        // _________________________________________
        // Factors

        //! Normalized Schwinger Electric field
        double norm_E_Schwinger;

        //! Inverse Normalized Schwinger Electric field
        double inv_norm_E_Schwinger;

        //! Espilon to check when tau is near 0
        const double epsilon_tau = 1e-100;

};

#endif
