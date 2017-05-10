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

#include "userFunctions.h"

// ---------------------------------------------------------------------------------------------------------------------
// Constructor for NLICompton
// ---------------------------------------------------------------------------------------------------------------------
NLICompton::NLICompton()
{
    Integfochi.resize(0);
}

// ---------------------------------------------------------------------------------------------------------------------
//
// Initialization of the parmeters for the nonlinear inverse Compton scattering
//
// 
// ---------------------------------------------------------------------------------------------------------------------
void NLICompton::initParams(Params& params)
{

    // If the namelist for Nonlinear Inverse Compton Scattering exists
    if( PyTools::nComponents("NLICompton") != 0 )
    {
        // Extraction of the parameter from the input file
        PyTools::extract("chipa_integfochi_min", chipa_integfochi_min, "NLICompton");
        PyTools::extract("chipa_integfochi_max", chipa_integfochi_max, "NLICompton");
        PyTools::extract("dim_integfochi", dim_integfochi, "NLICompton");
    }

    // Computation of the factor factor_dNphdt
    norm_lambda_compton = red_planck_cst*params.referenceAngularFrequency_SI/(electron_mass*c_vacuum);       
    factor_dNphdt = 1;

    // Some additional checks
    if (chipa_integfochi_min >= chipa_integfochi_max)
    {
        ERROR("chipa_integfochi_min (" << chipa_integfochi_min 
           << ") >= chipa_integfochi_max (" << chipa_integfochi_max << ")")
    } 
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
    i_chipa = int(floor(logchipa-chipa_integfochi_min)/delta_chipa_integfochi);

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
       logchipam = i_chipa*delta_chipa_integfochi + log10(chipa_integfochi_min);
       logchipap = logchipam + delta_chipa_integfochi;
   
       // Interpolation
       dNphdt = (Integfochi[i_chipa+1]*abs(logchipa-logchipam) + 
                 Integfochi[i_chipa]*abs(logchipap - logchipa))/delta_chipa_integfochi;
    }

    return factor_dNphdt*dNphdt*chipa/gfpa;

}

// ---------------------------------------------------------------------------------------------------------------------
// Computation of the table values of integfochi
// ---------------------------------------------------------------------------------------------------------------------
void NLICompton::compute_integfochi_table()
{
    double chipa; // Temporary particle chi value
    // For percentages
    int pct = 0;
    int dpct = 10;

    // Allocation of the array
    Integfochi.resize(dim_integfochi);

    // Computation of the delta
    delta_chipa_integfochi = (log10(chipa_integfochi_max) 
            - log10(chipa_integfochi_min))/(dim_integfochi-1);

    MESSAGE("Computation Integration F/chipa table:");

    // Loop over the table values
    for(unsigned int i = 0 ; i < dim_integfochi ; i++)
    {
        chipa = pow(i*delta_chipa_integfochi + log10(chipa_integfochi_min),10) ;

        Integfochi[i] = NLICompton::compute_integfochi(chipa,
                1e-40*chipa,0.98*chipa,400,1e-10);

        if (100.*i >= (double)(dim_integfochi*pct))
        {
            pct += dpct;
            MESSAGE(i + 1<< "/" << dim_integfochi << " - " << pct << "%");
        }
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
            << log10(chipa_integfochi_min) 
            << log10(chipa_integfochi_max);

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
//! using Gauss-Legendre for a given chipa value
//
//
//! \param chipa particle (electron for instance) quantum parameter
//! \param chipmin Minimal integration value (photon quantum parameter)
//! \param chipmax Maximal integration value (photon quantum parameter)
//! \param nbit number of points for integration
//! \param eps integration accuracy
// ---------------------------------------------------------------------------------------------------------------------
double NLICompton::compute_integfochi(double chipa,
        double chipmin,
        double chipmax,
        int nbit,
        double eps)
{

    //std::cout << "NLICompton::compute_integfochi" << std::endl;

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
        sync_emi = NLICompton::compute_sync_emissivity_ritus(chipa,chiph,200,1e-10);
        integ += gauleg_w[i]*sync_emi*log(10);
    }

    return integ;

}

// ---------------------------------------------------------------------------------------------------------------------
// Computation of the synchrotron emissivity following the formulae of Ritus
// ---------------------------------------------------------------------------------------------------------------------
double NLICompton::compute_sync_emissivity_ritus(double chipa, 
        double chiph, int nbit, double eps)
{

    //std::cout << "NLICompton::compute_sync_emissivity_ritus" << std::endl;

    // The photon quantum parameter should be below the electron one
    if (chipa > chiph)
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

        double y = chiph/(3.*chipa*(chipa-chiph));

        // Computation of Part. 1
        // Call the modified Bessel function to get K
        userFunctions::modified_bessel_IK(2./3.,2*y,I,dI,K,dK,50000,eps);

        part1 = (2. + 3.*chiph*y)*(K);

        // Computation of Part. 2
        // Using Gauss Legendre integration 

        userFunctions::gauss_legendre_coef(log10(2*y),log10(y)+50., gauleg_x, 
                gauleg_w, nbit, eps);

        part2 = 0;
        for(i=0 ; i< nbit ; i++)
        {
            y = pow(10.,gauleg_x[i]);
            userFunctions::modified_bessel_IK(1./3.,y,I,dI,K,dK,50000,eps);
            part2 += gauleg_w[i]*K*y*log(10.);
        }

        // Factor for final result
        y = 2*chiph/(3*pow(chipa,2.));

        return (part1 - part2)*y;


    }
    else if (chipa == chiph)
    {
        return 0;
    }
    else
    {
        ERROR("In compute_sync_emissivity_ritus: chipa " << chipa << " < chiph " << chiph);
    }    
}
