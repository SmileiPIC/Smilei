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

//  --------------------------------------------------------------------------------------------------------------------
//! Class Radiation
//  --------------------------------------------------------------------------------------------------------------------
class Radiation
{

    public:
        //! Creator for Radiation
        Radiation(Params& params, Species *species);
        virtual ~Radiation();

        //! Overloading of () operator
        virtual void operator() (Particles &particles,
                SmileiMPI* smpi,
                RadiationTables &RadiationTables,
                int istart,
                int iend,
                int ithread) = 0;

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
        double inline compute_chipa(double & charge_over_mass2,
                                     double & px, double & py, double & pz,
                                     double & gamma,
                                     double & Ex, double & Ey, double & Ez,
                                     double & Bx, double & By, double & Bz)
        {
            return fabs(charge_over_mass2)*inv_norm_E_Schwinger
                  * sqrt( fabs( pow(Ex*px + Ey*py + Ez*pz,2)
                  - pow(gamma*Ex - By*pz + Bz*py,2)
                  - pow(gamma*Ey - Bz*px + Bx*pz,2)
                  - pow(gamma*Ez - Bx*py + By*px,2)));
        };

        //! Return the total normalized radiated energy
        double inline getRadiatedEnergy()
        {
            return radiated_energy;
        };

        //! Computation of the quantum parameter for the given
        //! thread of particles
        //! \param Particles class containg the particle property arrays
        //! \param smpi class for mpi parameters
        //! \param istart      Index of the first particle
        //! \param iend        Index of the last particle
        //! \param ithread     Thread index
        void compute_thread_chipa(Particles &particles,
                SmileiMPI* smpi,
                int istart,
                int iend,
                int ithread);


    protected:

        // ________________________________________
        // General parameters

        //! Dimension of position
        int nDim_;

        //! Inverse species mass
        double one_over_mass_;

        //! Time step
        double dt;

        //! Radiated energy of the total thread
        double radiated_energy;

        // _________________________________________
        // Factors

        //! Fine structure constant
        const double fine_struct_cst = 7.2973525698e-3;

        //! Reduced Planck Constant (J.s)
        const double red_planck_cst = 1.054571628E-34;

        //! Electron mass (kg)
        const double electron_mass = 9.109382616e-31;

        //! Speed of light in vacuum (m/s)
        const double c_vacuum = 299792458;

        //! Normalized Schwinger Electric field
        double norm_E_Schwinger;

        //! Inverse Normalized Schwinger Electric field
        double inv_norm_E_Schwinger;

    private:

};//END class

#endif
