#ifndef BINARYPROCESSDATA_H
#define BINARYPROCESSDATA_H

#ifdef SMILEI_ACCELERATOR_GPU
#define SMILEI_BINARYPROCESS_BUFFERSIZE 32
#else
#define SMILEI_BINARYPROCESS_BUFFERSIZE 8
#endif

#define SMILEI_BINARYPROCESS_FLOAT float

#include "Particles.h"

//! Contains the relativistic kinematic quantities associated to the collision of two particles noted 1 and 2
struct BinaryProcessData
{
    
    //! Indices of both particles in their Particles object
    uint32_t i[2][SMILEI_BINARYPROCESS_BUFFERSIZE];
    
    //! Particles objects of both particles
    Particles * p[2][SMILEI_BINARYPROCESS_BUFFERSIZE];
    
    //! Species index
    uint32_t ispec[2][SMILEI_BINARYPROCESS_BUFFERSIZE];
    
    //! Masses and mass ratio
    SMILEI_BINARYPROCESS_FLOAT m[2][SMILEI_BINARYPROCESS_BUFFERSIZE], R[SMILEI_BINARYPROCESS_BUFFERSIZE];
    
    //! Weights
    double W[2][SMILEI_BINARYPROCESS_BUFFERSIZE];

    //! Charges
    short q[2][SMILEI_BINARYPROCESS_BUFFERSIZE];
    
    //! Momenta
    double px[2][SMILEI_BINARYPROCESS_BUFFERSIZE], py[2][SMILEI_BINARYPROCESS_BUFFERSIZE], pz[2][SMILEI_BINARYPROCESS_BUFFERSIZE];
    
    //! Correction to apply to the cross-sections due to the difference in weight
    SMILEI_BINARYPROCESS_FLOAT dt_correction[SMILEI_BINARYPROCESS_BUFFERSIZE];
    
    //! Sum of both momenta
    SMILEI_BINARYPROCESS_FLOAT px_tot[SMILEI_BINARYPROCESS_BUFFERSIZE], py_tot[SMILEI_BINARYPROCESS_BUFFERSIZE], pz_tot[SMILEI_BINARYPROCESS_BUFFERSIZE];
    
    //! Lorentz invariant = energy of one particle in the frame of the other
    SMILEI_BINARYPROCESS_FLOAT gamma0[SMILEI_BINARYPROCESS_BUFFERSIZE];
    
    //! Momentum of the particles expressed in the COM frame
    SMILEI_BINARYPROCESS_FLOAT px_COM[SMILEI_BINARYPROCESS_BUFFERSIZE], py_COM[SMILEI_BINARYPROCESS_BUFFERSIZE], pz_COM[SMILEI_BINARYPROCESS_BUFFERSIZE];
    SMILEI_BINARYPROCESS_FLOAT p_COM[SMILEI_BINARYPROCESS_BUFFERSIZE];
    
    //! Lorentz factors
    SMILEI_BINARYPROCESS_FLOAT gamma[2][SMILEI_BINARYPROCESS_BUFFERSIZE], gamma_tot[SMILEI_BINARYPROCESS_BUFFERSIZE];
    //! Lorentz factors expressed in the COM frame
    SMILEI_BINARYPROCESS_FLOAT gamma_COM0[SMILEI_BINARYPROCESS_BUFFERSIZE], gamma_tot_COM[SMILEI_BINARYPROCESS_BUFFERSIZE];
    
    //! Relative velocity
    SMILEI_BINARYPROCESS_FLOAT vrel[SMILEI_BINARYPROCESS_BUFFERSIZE];
    
    //! Whether the first species is electron
    bool electronFirst;
    
    //! Debye length
    SMILEI_BINARYPROCESS_FLOAT debye;
    
    //! Thomas-Fermi length
    SMILEI_BINARYPROCESS_FLOAT lTF[SMILEI_BINARYPROCESS_BUFFERSIZE];
    
    //! Product of atomic numbers (for e-i screening)
    SMILEI_BINARYPROCESS_FLOAT Z1Z2[SMILEI_BINARYPROCESS_BUFFERSIZE];
    
    //! Which group of species may be the screened species
    uint32_t screening_group;
    
    SMILEI_BINARYPROCESS_FLOAT n123, n223;
    
    //! Other buffers available for the binary processes
    SMILEI_BINARYPROCESS_FLOAT buffer1[SMILEI_BINARYPROCESS_BUFFERSIZE];
    SMILEI_BINARYPROCESS_FLOAT buffer2[SMILEI_BINARYPROCESS_BUFFERSIZE];
    SMILEI_BINARYPROCESS_FLOAT buffer3[SMILEI_BINARYPROCESS_BUFFERSIZE];
    SMILEI_BINARYPROCESS_FLOAT buffer4[SMILEI_BINARYPROCESS_BUFFERSIZE];
    SMILEI_BINARYPROCESS_FLOAT buffer5[SMILEI_BINARYPROCESS_BUFFERSIZE];
    
};

#endif
