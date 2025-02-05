// ----------------------------------------------------------------------------
//! \file MultiphotonBreitWheeler.h
//
//! \brief This file contains the class functions for the generic class
//!  MultiphotonBreitWheeler for the photon decay into pairs via the multiphoton
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
#include "Params.h"
#include "Random.h"

//  ----------------------------------------------------------------------------
//! Class Radiation
//  ----------------------------------------------------------------------------
class MultiphotonBreitWheeler
{
public:

    //! Creator for Radiation
    MultiphotonBreitWheeler( Params &params, Species *species, Random * rand );
    ~MultiphotonBreitWheeler();

    //! Overloading of () operator
    //! \param particles   particle object containing the particle
    //!                    properties of the current species
    //! \param smpi        MPI properties
    //! \param MultiphotonBreitWheelerTables Cross-section data tables and useful functions
    //!                     for multiphoton Breit-Wheeler
    //! \param pair_energy energy converted into pairs
    //! \param istart      Index of the first particle
    //! \param iend        Index of the last particle
    //! \param ithread     Thread index
    void operator()( Particles &particles,
                     SmileiMPI* smpi,
                     Particles** new_pair,
                     Species ** new_pair_species,
                     MultiphotonBreitWheelerTables &mBW_tables,
                     double & pair_energy,
                     int istart,
                     int iend,
                     int ithread, int ibin = 0, int ipart_ref = 0 );
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
    inline double __attribute__((always_inline)) computePhotonChi(
                                 double kx, double ky, double kz,
                                 double gamma,
                                 double Ex, double Ey, double Ez,
                                 double Bx, double By, double Bz )
    {

        return inv_norm_E_Schwinger_
               * std::sqrt( std::fabs( ( Ex*kx + Ey*ky + Ez*kz) * (Ex*kx + Ey*ky + Ez*kz)
                             - ( gamma*Ex - By*kz + Bz*ky ) * ( gamma*Ex - By*kz + Bz*ky )
                             - ( gamma*Ey - Bz*kx + Bx*kz ) * ( gamma*Ey - Bz*kx + Bx*kz )
                             - ( gamma*Ez - Bx*ky + By*kx ) * ( gamma*Ez - Bx*ky + By*kx ) ) );
    };

    //! Computation of the quantum parameter for the given
    //! thread of photons
    //! \param Particles class containg the particle property arrays
    //! \param smpi class for mpi parameters
    //! \param istart      Index of the first particle
    //! \param iend        Index of the last particle
    //! \param ithread     Thread index
    void computeThreadPhotonChi( Particles &particles,
                               SmileiMPI *smpi,
                               int istart,
                               int iend,
                               int ithread, int ipart_ref = 0 );

    //! Clean photons that decayed into pairs (weight <= 0)
    //! \param particles   particle object containing the particle
    //!                    properties of the current species
    //! \param smpi        MPI properties
    //! \param ibin        Index of the current bin
    //! \param nbin        Number of bins
    //! \param bmin        Pointer toward the first particle index of the bin in the Particles object
    //! \param bmax        Pointer toward the last particle index of the bin in the Particles object
    //! \param ithread     Thread index
    void removeDecayedPhotons(
        Particles &particles,
        SmileiMPI *smpi,
        int ibin, int nbin,
        int *bmin, int *bmax, int ithread );


    //! Clean photons that decayed into pairs (weight <= 0) and resize each bin
    //! But keeping the space between bins (so called no compression)
    //! Developers have to be aware that the space exists using the Particles bin indexes
    //! \param particles   particle object containing the particle
    //!                    properties of the current species
    //! \param smpi        MPI properties
    //! \param ibin        Index of the current bin
    //! \param bmin        Pointer toward the first particle index of the bin in the Particles object
    //! \param bmax        Pointer toward the last particle index of the bin in the Particles object
    //! \param ithread     Thread index
//#ifdef SMILEI_ACCELERATOR_GPU_OACC
//    #pragma acc routine seq
//#endif
    void removeDecayedPhotonsWithoutBinCompression(
        Particles &particles,
        SmileiMPI *smpi,
        int ibin,
        int *bmin, int *bmax, int ithread );

    //! Return the sampling for each pair
    int getPairCreationSampling(int i) {
        return mBW_pair_creation_sampling_[i];
    }

    //! Return the pair converted energy
    // double inline getPairEnergy( void )
    // {
    //     return pair_converted_energy_;
    // }

    // Local array of new pairs of electron-positron
    // Particles new_pair[2];

    // Local array of new pairs of electron-positron per bin
    std::vector<Particles *> new_pair_per_bin;


private:

    // ________________________________________
    // General parameters

    //! Dimension of position
    int n_dimensions_;

    //! Time step
    double dt_;

    // Number of pairs created per event
    int mBW_pair_creation_sampling_[2];

    // Inverse of the number of pairs created per even
    double mBW_pair_creation_inv_sampling_[2];

    //! Threshold under which pair creation is not considered
    double chiph_threshold_;

    //! Local random generator
    Random * rand_;

    // _________________________________________
    // Factors

    //! Normalized Schwinger Electric field
    double norm_E_Schwinger_;

    //! Inverse Normalized Schwinger Electric field
    double inv_norm_E_Schwinger_;

    //! Espilon to check when tau is near 0
    static constexpr double epsilon_tau_ = 1e-100;

};

#endif
