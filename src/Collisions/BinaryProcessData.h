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
    double m1, m2, m12;
    
    //! Minimum / maximum weight
    double minW, maxW;
    
    //! Whether the first species is electron
    bool electronFirst;
    
    //! Correction to apply to the cross-sections due to the difference in weight
    double dt_correction;
    
    //! Velocity of the Center-Of-Mass, expressed in the lab frame
    double COM_vx, COM_vy, COM_vz;
    
    //! Lorentz factor of the COM, expressed in the lab frame
    double COM_gamma;
    
    //! Momentum of the particles expressed in the COM frame
    double px_COM, py_COM, pz_COM, p_COM;
    
    //! Lorentz factors
    double gamma1, gamma2;
    //! Lorentz factors expressed in the COM frame
    double gamma1_COM, gamma2_COM;
    
    //! Relative velocity
    double vrel, vrel_corr;
    
    //! Debye length squared
    double debye2;
    
    double term1, term3, term5, n123, n223;
};

#endif
