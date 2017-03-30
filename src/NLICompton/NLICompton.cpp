/********************************************************************

 CPP Nonlinear Inverse Compton Scattering

 This header contains the definition of the class NLICompton.

********************************************************************/

#include "NLICompton.h" 

#include <cstring>
#include <iostream>

#include <cmath>

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
