// ----------------------------------------------------------------------------
//! \file NLICompton.cpp
//
//! \brief Nonlinear Inverse Compton Scattering
//
//! \details This header contains the definition of the class NLICompton.
//
// ----------------------------------------------------------------------------

#include "NLICompton.h" 

#include <cstring>
#include <iostream>
#include <fstream>

#include <cmath>

#include "Params.h"
#include "userFunctions.h"

// ---------------------------------------------------------------------------------------------------------------------
// Constructor for NLICompton
// ---------------------------------------------------------------------------------------------------------------------
NLICompton::NLICompton()
{
    Integfochi.resize(0);
}

// ---------------------------------------------------------------------------------------------------------------------
// Computation of the Cross Section dNph/dt which is also 
// the number of photons generated per time unit.
// ---------------------------------------------------------------------------------------------------------------------
double NLICompton::compute_dNphdt(double chipa,double gfpa)
{

    // Log of the particle quantum parameter chipa
    double logchipa;
    double logchipam, logchipap;
    // Index
    unsigned int i_chipa;
    // final value
    double dNphdt;

    logchipa = log10(chipa);

    // Lower index for interpolation in the table integfochi
    i_chipa = int(floor(logchipa-chie_integfochi_min)/delta_chie_integfochi);

    // If we are not in the table...
    if (i_chipa < 0)
    {
        i_chipa = 0;
        dNphdt = Integfochi[i_chipa];
    }
    else if (i_chipa >= dim_integfochi-1)
    {
        i_chipa = dim_integfochi-2;
        dNphdt = Integfochi[i_chipa];
    }
    else
    {
       // Upper and minor values for linear interpolation
       logchipam = i_chipa*delta_chie_integfochi + log10(chie_integfochi_min);
       logchipap = logchipam + delta_chie_integfochi;
   
       // Interpolation
       dNphdt = (Integfochi[i_chipa+1]*abs(logchipa-logchipam) + 
                 Integfochi[i_chipa]*abs(logchipap - logchipa))/delta_chie_integfochi;
    }

    return factor_dNphdt*dNphdt*chipa/gfpa;

}

// ---------------------------------------------------------------------------------------------------------------------
// Computation of the table values of integfochi
// ---------------------------------------------------------------------------------------------------------------------
void NLICompton::compute_integfochi_table()
{
    double chie; // Temporary particle chi value

    // Allocation of the array
    Integfochi.resize(dim_integfochi);

    // Computation of the delta
    delta_chie_integfochi = (log10(chie_integfochi_max) 
            - log10(chie_integfochi_min))/(dim_integfochi-1);

    // Loop over the table values
    for(unsigned int i = 0 ; i < dim_integfochi ; i++)
    {
        chie = pow(i*delta_chie_integfochi + log10(chie_integfochi_min),10) ;

        Integfochi[i] = NLICompton::compute_integfochi(chie,
                1e-40*chie,0.99*chie,400,1e-15);

    }  

}


// ---------------------------------------------------------------------------------------------------------------------
// Ouput in a file of the table values of integfochi
// ---------------------------------------------------------------------------------------------------------------------
void NLICompton::output_integfochi_table()
{

    std::ofstream file;
    file.open("tab_integfochi.bin",std::ios::binary);

    if (file.is_open()) {

        file << dim_integfochi 
            << log10(chie_integfochi_min) 
            << log10(chie_integfochi_max);

        // Loop over the table values
        for(unsigned int i = 0 ; i < dim_integfochi ; i++)
        {
            file << Integfochi[i];
        }

    file.close();
    }
}


// ---------------------------------------------------------------------------------------------------------------------
//! \brief 
//! Compute integration of F/chi between 
//! using Gauss-Legendre for a given chie value
//
//
//! \param chie particle (electron for instance) quantum parameter
//! \param chipmin Minimal integration value (photon quantum parameter)
//! \param chipmax Maximal integration value (photon quantum parameter)
//! \param nbit number of points for integration
//! \param eps integration accuracy
// ---------------------------------------------------------------------------------------------------------------------
double NLICompton::compute_integfochi(double chie,
        double chipmin,
        double chipmax,
        int nbit,
        double eps)
{

    // Arrays for Gauss-Legendre integration
    double * gauleg_x = new double[nbit];
    double * gauleg_w = new double[nbit];
    // Photon quantum parameter
    double chiph;
    // Integration result
    double integ;
    // Synchrotron emissivity
    double sync_emi;
    // Iterator
    int i;

    // gauss Legendre coefficients
    userFunctions::gauss_legendre_coef(log10(chipmin),log10(chipmax), 
            gauleg_x, gauleg_w, nbit, eps);

    // Integration loop
    integ = 0;
    for(i=0 ; i< nbit ; i++)
    {
        chiph = pow(10.,gauleg_x[i]);
        sync_emi = NLICompton::compute_sync_emissivity_ritus(chie,chiph,200,1e-15);
        integ += gauleg_w[i]*sync_emi*log(10);
    }

    return integ;

}

// ---------------------------------------------------------------------------------------------------------------------
// Computation of the synchrotron emissivity following the formulae of Ritus
// ---------------------------------------------------------------------------------------------------------------------
double NLICompton::compute_sync_emissivity_ritus(double chie, 
        double chiph, int nbit, double eps)
{
    // The photon quantum parameter should be below the electron one
    if (chie > chiph)
    {
        // Arrays for Gauss-Legendre integration
        double * gauleg_x = new double[nbit];
        double * gauleg_w = new double[nbit];
        // Values for Bessel results       
        double I,dI;
        double K,dK;
        // Parts of the formulae
        double part1,part2;
        // Iterator
        int i;

        double y = chiph/(3.*chie*(chie-chiph));

        // Computation of Part. 1
        // Call the modified Bessel function to get K
        userFunctions::modified_bessel_IK(2./3.,y,I,dI,K,dK,20000,eps);

        part1 = (2. + 3.*chiph*y)*(K);

        // Computation of Part. 2
        // Using Gauss Legendre integration 

        userFunctions::gauss_legendre_coef(log10(2*y),log10(y)+50., gauleg_x, 
                gauleg_w, nbit, eps);

        part2 = 0;
        for(i=0 ; i< nbit ; i++)
        {
            y = pow(10.,gauleg_x[i]);
            userFunctions::modified_bessel_IK(1./3.,y,I,dI,K,dK,20000,eps);
            part2 += gauleg_w[i]*K*y*log(10.);
        }

        // Factor for final result
        y = 2*chiph/(3*pow(chie,2.));

        return (part1 - part2)*y;


    }
    else if (chie == chiph)
    {
        return 0;
    }
    else
    {
        ERROR("In compute_sync_emissivity_ritus: chie " << chie << " < chiph " << chiph);
    }    
}
