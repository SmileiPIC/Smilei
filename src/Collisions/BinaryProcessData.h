#ifndef BINARYPROCESSDATA_H
#define BINARYPROCESSDATA_H

#include "Particles.h"

//! Contains the relativistic kinematic quantities associated to the collision of two particles noted 1 and 2
struct BinaryProcessData
{
    //! Number of particles in buffer
    size_t n;
    
    //! Indices of both particles in their Particles object
    std::vector<size_t> i[2];
    
    //! Particles objects of both particles
    std::vector<Particles *> p[2];
    
    //! Masses
    std::vector<double> m[2];
    
    //! Weights
    std::vector<double> W[2];

    //! Charges
    std::vector<short> q[2];
    
    //! Correction to apply to the cross-sections due to the difference in weight
    std::vector<double> dt_correction;
    
    //! Momenta
    std::vector<double> px[2], py[2], pz[2];
    
    //! Sum of both momenta
    std::vector<double> px_tot, py_tot, pz_tot;
    
    //! Lorentz invariant = energy of one particle in the frame of the other
    std::vector<double> gamma0;
    
    //! Momentum of the particles expressed in the COM frame
    std::vector<double> px_COM, py_COM, pz_COM;
    std::vector<double> p_COM;
    
    //! Lorentz factors
    std::vector<double> gamma[2], gamma_tot;
    //! Lorentz factors expressed in the COM frame
    std::vector<double> gamma_COM0, gamma_tot_COM;
    
    //! Relative velocity
    std::vector<double> vrel, vrel_corr;
    
    //! Whether the first species is electron
    bool electronFirst;
    
    //! Debye length squared
    double debye;
    
    //! Thomas-Fermi length
    std::vector<double> lTF;
    
    //! Product of atomic numbers (for e-i screening)
    std::vector<double> Z1Z2;
    
    double n123, n223;
    
    void resize( size_t N ) {
        i[0].resize( N ); i[1].resize( N );
        p[0].resize( N ); p[1].resize( N );
        m[0].resize( N ); m[1].resize( N );
        W[0].resize( N ); W[1].resize( N );
        q[0].resize( N ); q[1].resize( N );
        px[0].resize( N ); py[0].resize( N ); pz[0].resize( N );
        px[1].resize( N ); py[1].resize( N ); pz[1].resize( N );
        dt_correction.resize( N );
        px_tot.resize( N ); py_tot.resize( N ); pz_tot.resize( N );
        gamma0.resize( N );
        px_COM.resize( N ); py_COM.resize( N ); pz_COM.resize( N );
        p_COM.resize( N );
        gamma[0].resize( N ); gamma[1].resize( N ); gamma_tot.resize( N );
        gamma_COM0.resize( N ); gamma_tot_COM.resize( N );
        vrel.resize( N ); vrel_corr.resize( N );
        lTF.resize( N );
        Z1Z2.resize( N );
    }
};

#endif
