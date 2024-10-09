#ifndef BINARYPROCESSDATA_H
#define BINARYPROCESSDATA_H

#include "Particles.h"

//! Contains the relativistic kinematic quantities associated to the collision of two particles noted 1 and 2
struct BinaryProcessData
{
    //! Particles objects for both macro-particles
    Particles *p1, *p2;
    
    //! Indices of both particles
    unsigned int i1, i2;
    
    //! Masses
    double m1, m2, m21;
    
    //! Minimum / maximum weight
    double minW, maxW;
    
    //! Whether the first species is electron
    bool electronFirst;
    
    //! Correction to apply to the cross-sections due to the difference in weight
    double dt_correction;
    
    //! Sum of both momenta
    double px_tot, py_tot, pz_tot;
    
    //! Lorentz invariant = energy of one particle in the frame of the other
    double E0;
    
    //! Momentum of the particles expressed in the COM frame
    double px_COM, py_COM, pz_COM;
    double p2_COM, p_COM;
    
    //! Lorentz factors
    double gamma1, gamma2, gamma_tot;
    //! Lorentz factors expressed in the COM frame
    double gamma1_COM, gamma2_COM, gamma_tot_COM;
    
    //! Relative velocity
    double vrel, vrel_corr, p_gamma_COM;
    
    //! Debye length squared
    double debye;
    
    //! Thomas-Fermi length
    double lTF;
    
    //! Product of atomic numbers (for e-i screening)
    double Z1Z2;
    
    double n123, n223;
};

#endif
