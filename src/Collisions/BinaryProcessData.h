#ifndef BINARYPROCESSDATA_H
#define BINARYPROCESSDATA_H

#ifdef SMILEI_ACCELERATOR_GPU
#define SMILEI_BINARYPROCESS_BUFFERSIZE 32
#else
#define SMILEI_BINARYPROCESS_BUFFERSIZE 8
#endif

#include "Particles.h"

//! Contains the relativistic kinematic quantities associated to the collision of two particles noted 1 and 2
struct BinaryProcessData
{
    
    //! Indices of both particles in their Particles object
    size_t i[2][SMILEI_BINARYPROCESS_BUFFERSIZE];
    
    //! Particles objects of both particles
    Particles * p[2][SMILEI_BINARYPROCESS_BUFFERSIZE];
    
    //! Species index
    size_t ispec[2][SMILEI_BINARYPROCESS_BUFFERSIZE];
    
    //! Masses and mass ratio
    double m[2][SMILEI_BINARYPROCESS_BUFFERSIZE], R[SMILEI_BINARYPROCESS_BUFFERSIZE];
    
    //! Weights
    double W[2][SMILEI_BINARYPROCESS_BUFFERSIZE];

    //! Charges
    short q[2][SMILEI_BINARYPROCESS_BUFFERSIZE];
    
    //! Correction to apply to the cross-sections due to the difference in weight
    double dt_correction[SMILEI_BINARYPROCESS_BUFFERSIZE];
    
    //! Momenta
    double px[2][SMILEI_BINARYPROCESS_BUFFERSIZE], py[2][SMILEI_BINARYPROCESS_BUFFERSIZE], pz[2][SMILEI_BINARYPROCESS_BUFFERSIZE];
    
    //! Sum of both momenta
    double px_tot[SMILEI_BINARYPROCESS_BUFFERSIZE], py_tot[SMILEI_BINARYPROCESS_BUFFERSIZE], pz_tot[SMILEI_BINARYPROCESS_BUFFERSIZE];
    
    //! Lorentz invariant = energy of one particle in the frame of the other
    double gamma0[SMILEI_BINARYPROCESS_BUFFERSIZE];
    
    //! Momentum of the particles expressed in the COM frame
    double px_COM[SMILEI_BINARYPROCESS_BUFFERSIZE], py_COM[SMILEI_BINARYPROCESS_BUFFERSIZE], pz_COM[SMILEI_BINARYPROCESS_BUFFERSIZE];
    double p_COM[SMILEI_BINARYPROCESS_BUFFERSIZE];
    
    //! Lorentz factors
    double gamma[2][SMILEI_BINARYPROCESS_BUFFERSIZE], gamma_tot[SMILEI_BINARYPROCESS_BUFFERSIZE];
    //! Lorentz factors expressed in the COM frame
    double gamma_COM0[SMILEI_BINARYPROCESS_BUFFERSIZE], gamma_tot_COM[SMILEI_BINARYPROCESS_BUFFERSIZE];
    
    //! Relative velocity
    double vrel[SMILEI_BINARYPROCESS_BUFFERSIZE];
    
    //! Whether the first species is electron
    bool electronFirst;
    
    //! Debye length squared
    double debye;
    
    //! Thomas-Fermi length
    double lTF[SMILEI_BINARYPROCESS_BUFFERSIZE];
    
    //! Product of atomic numbers (for e-i screening)
    double Z1Z2[SMILEI_BINARYPROCESS_BUFFERSIZE];
    
    //! Which group of species may be the screened species
    int screening_group;
    
    double n123, n223;
    
    //! Other buffers available for the binary processes
    double buffer1[SMILEI_BINARYPROCESS_BUFFERSIZE];
    double buffer2[SMILEI_BINARYPROCESS_BUFFERSIZE];
    double buffer3[SMILEI_BINARYPROCESS_BUFFERSIZE];
    double buffer4[SMILEI_BINARYPROCESS_BUFFERSIZE];
    double buffer5[SMILEI_BINARYPROCESS_BUFFERSIZE];
    
};

#endif
