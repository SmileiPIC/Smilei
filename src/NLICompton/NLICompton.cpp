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
// Computation of the table values of integfochi
// ---------------------------------------------------------------------------------------------------------------------
void NLICompton::compute_integfochi()
{
    double chie; // Temporary chi value

    // Allocation of the array
    Integfochi.resize(dim_integfochi);

    // Computation of the delta
    delta_chie_integfochi = (log10(chie_integfochi_max) - log10(chie_integfochi_min))/(dim_integfochi-1);

    // Loop over the table values
    for(unsigned int i = 0 ; i < dim_integfochi ; i++)
    {
        chie = pow(i*delta_chie_integfochi + log10(chie_integfochi_min),10) ;
    }  

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

       userFunctions::gauss_legendre_coef(log(2*y),log(y)+50., gauleg_x, 
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
