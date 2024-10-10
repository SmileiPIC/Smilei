#ifndef BINARYPROCESSDATA_H
#define BINARYPROCESSDATA_H

#include "Particles.h"

//! Contains the relativistic kinematic quantities associated to the collision of two particles noted 1 and 2
struct BinaryProcessData
{
    static constexpr size_t max_buffer_size_ = 8;
    
    //! Number of particles in buffer
    size_t n;
    
    //! Indices of both particles in their Particles object
    size_t i[2][max_buffer_size_];
    
    //! Particles objects of both particles
    Particles * p[2][max_buffer_size_];
    
    //! Masses
    double m[2][max_buffer_size_];
    
    //! Weights
    double W[2][max_buffer_size_];

    //! Charges
    short q[2][max_buffer_size_];
    
    //! Correction to apply to the cross-sections due to the difference in weight
    double dt_correction[max_buffer_size_];
    
    //! Momenta
    double px[2][max_buffer_size_], py[2][max_buffer_size_], pz[2][max_buffer_size_];
    
    //! Sum of both momenta
    double px_tot[max_buffer_size_], py_tot[max_buffer_size_], pz_tot[max_buffer_size_];
    
    //! Lorentz invariant = energy of one particle in the frame of the other
    double gamma0[max_buffer_size_];
    
    //! Momentum of the particles expressed in the COM frame
    double px_COM[max_buffer_size_], py_COM[max_buffer_size_], pz_COM[max_buffer_size_];
    double p_COM[max_buffer_size_];
    
    //! Lorentz factors
    double gamma[2][max_buffer_size_], gamma_tot[max_buffer_size_];
    //! Lorentz factors expressed in the COM frame
    double gamma_COM0[max_buffer_size_], gamma_tot_COM[max_buffer_size_];
    
    //! Relative velocity
    double vrel[max_buffer_size_], vrel_corr[max_buffer_size_];
    
    //! Whether the first species is electron
    bool electronFirst;
    
    //! Debye length squared
    double debye;
    
    //! Thomas-Fermi length
    double lTF[max_buffer_size_];
    
    //! Product of atomic numbers (for e-i screening)
    double Z1Z2[max_buffer_size_];
    
    double n123, n223;
};

#endif
