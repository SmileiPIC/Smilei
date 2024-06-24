#include "continuity_tool.h" // declares double au_to_eV;
#include "Ionization.h"
#include "IonizationTunnelBSI.h"
#include "IonizationTables.h" // can use the raw ionizationEnergy[100/2*101] (all in eV), but also the function ionization_energy(int atomic_number, int Zstar) which returns: ionizationEnergy[( atomic_number*(atomic_number-1) )/2 + Zstar]
// just write IonizationTables::ionizationEnergy[atomic_number_*(atomic_number-1)/2 + Zstar] to get IonizPotential (in eV) of Element with atmomic # = atomic_number and charge state Zstar, i.e. Element^{+Zstar}
// For instance IonizationTables::ionizationEnergy[13*12/2 + 3] is the ionization energy of Al3+ in eV.

#include <cmath>
#include <vector>

double au_to_eV = 27.2116; // will be used by continuity_tool( , , ) to convert Potential[Z] (in AtomicUnits) to eV.


int continuity_tool(unsigned int Zp1, double E, double alpha, double beta, double gamma, double E_cr, double Potential, unsigned int atomic_number_, double au_to_w0) {
    // returns:
    // 0 for when it chooses Tunnel ADK-PPT Ionization Rate, 
    // 1 for when it chooses the BSI-Quadratic Ionization Rate, 
    // 2 for when it chooses BSI-Linear Ionization Rate

    // Aim is to implement the 3-piece general rate, so a method to choose between Tunnel or BSI linear or BSI quadratic.
    // E_cr from the arguments of this function is useless at the moment. 
    // alpha, beta, gamma coefficients are single numbers, for ion with charge state Z = Zp1 - 1, for a fixed atomic number (of the species under consideration). Z=0 is a neutral atom, Z=1 is singly ionized ion, i.e. Na^{+1}.
    // alpha = alpha_tunnel[Z], beta = beta_tunnel[Z], gamma = gamma_tunnel[Z], where Z = Zp1 - 1 and 
    // where alpha_tunnel[Z] uses Z and Potential[Z] (IonizationPotential (in AtomicUnits) of Element with atomic number = atomic_number_ and charge state Z (Z=0 means neutral atom))
    // (double) Potential from arguments is the IP in AtomicUnits for current charge state Z (Z=0 is neutral atom) for a fixed atomic no = atomic_number_.
    // Potential from the arguments is not actually needed!!! 
    unsigned int Z = Zp1 - 1;

    double Potential_in_eV = Potential * au_to_eV; // AtomicUnits to eV conversion on Potential from arguments of this function.  // not actually needed !!
    // Potential_in_eV will be used for BSI linear rate implemented below. There we need the ratio of Potential_in_eV[] and 13.6 eV  
    // not actually needed anymore!
    
    double BSI_rate_quadratic = BSI_rate_quadratic_function(Zp1, E, atomic_number_, au_to_w0); // returns rate for charge state Z and atomic_number=atomic_number_
    double BSI_rate_linear = BSI_rate_linear_function(Zp1, E, atomic_number_, au_to_w0); // returns rate for charge state Z and atomic_number=atomic_number_
    double Tunnel_rate = Tunnel_rate_function(Z, E, alpha, beta, gamma); // Z = Zp1 - 1;

    //double BSI_rate_quadratic_loglogder_wrtE = BSI_rate_quadratic_function_loglogder_wrtE(Zp1, E, atomic_number_);
    //double BSI_rate_linear_loglogder_wrtE = BSI_rate_linear_function_loglogder_wrtE(Zp1, E, atomic_number_);
    //double Tunnel_rate_loglogder_wrtE = Tunnel_rate_function_loglogder_wrtE(Z, E, alpha, beta, gamma);

// Tony Arber 2015, the review paper on PIC codes. Arber_2015_Plasma_Phys._Control._Fusion_57_113001
// There he explains how to make the PIC code choose, at each timestep, between Tunnel and BSI.
// So he tells us how to implement a 2-piece formula for the ionization rate.

//     if (BSI_rate_quadratic_loglogder_wrtE > Tunnel_rate_loglogder_wrtE) { // Conditions for a 2-piece general rate, where the aim is to choose between Tunnel and BSI quadratic rate.
//         return true; // i.e. use BSI rate because E > Eprime (location of the intersection of the 2 rates)
//     } 
//     else if (BSI_rate_quadratic > Tunnel_rate) { // simulates: else return minimum(BSI_rate_quadratic, Tunnel_rate)
//         return false; // i.e. use Tunnel Rate 
//     }    
//     else return true; // i.e. use BSI Rate.

// Below I implement a method to choice from the 3 rates (i.e. the 3 piece formula)
// based on something simpler than in Arber where he uses log-log derivatives.
// Just look at the graphs from Golovanov 2020 paper (arxhiv one for nicely coloured graphs).
    if (BSI_rate_quadratic >= BSI_rate_linear) { // all these rates which are compared here are in SMILEI units
        return 2; // 2 means return BSI-Linear rate
    }
    else if (std::min(Tunnel_rate, BSI_rate_quadratic) == BSI_rate_quadratic) { // this returns min(Tunnel, BSI-Quadratic)
        return 1; // 1 means return BSI-Quadratic rate
    }
    else return 0; // 0 means return Tunnel rate. The min is Tunnel rate
    
}

   
// TUNNEL Ionization Rates
// #############################################################
// 1) ADK-PPT static rate in units of atomic_unit_of_freq, all physical quantities inside rate formula are in atomic units.
double Tunnel_rate_function(unsigned int Z, double E, double alpha, double beta, double gamma) { // E is in Atomic Units
    // Z from arguments is useless here as alpha, beta, gamma are for that Z already.
    double one_third = 1.0/3.0;
    double invE = 1./E ;
    double delta = gamma * invE;
    
    double Tunnel_rate = beta * exp( -delta*one_third + alpha*log(delta) ); // E depenedency comes through delta
    return Tunnel_rate; // returned rate is in SMILEI UNITS due to beta from above
}

// double Tunnel_rate_function_der_wrtE(unsigned int Z, double E, double alpha, double beta, double gamma) {
//     // I wrote the expression for derivative using gamma only, nothing related to delta.
//     double one_third = 1.0/3.0;
//     double invE = 1./E ;
//     // double delta = gamma * invE; // I wrote the expression for derivative using gamma only.
    
//     double Tunnel_rate_der_wrtE = beta * exp(- gamma*one_third*invE) * pow((gamma*invE), alpha) * (gamma*one_third*pow(invE, 2) - alpha*invE);
//     return Tunnel_rate_der_wrtE;
// }

// double Tunnel_rate_function_loglogder_wrtE(unsigned int Z, double E, double alpha, double beta, double gamma) {
//     double Tunnel_rate = Tunnel_rate_function(Z, E, alpha, beta, gamma);
//     double Tunnel_rate_der_wrtE = Tunnel_rate_function_der_wrtE(Z, E, alpha, beta, gamma);

//     double Tunnel_rate_loglogder_wrtE = (E/Tunnel_rate) * Tunnel_rate_der_wrtE; // log-log derivative is: (x/y) * (dy/dx) where in our case: y = Rate, x = E-field.
//     return Tunnel_rate_loglogder_wrtE;
// }

// BSI rates
// ##############################################################################
// BSI rates below are calculated in atomic units, for both the Quadratic rate and the Linear rate.
// 1) QUADRATIC dependence on local E-field <--> Bauer and Mulser, in PhysRevA Vol.59, No.1, 1999.
double BSI_rate_quadratic_function(unsigned int Zp1, double E, unsigned int atomic_number_, double au_to_w0) { // E-field E from arguments has to be in AtomicUnits.
    // double BSI_rate_quadratic = 2.4 * (pow(E,2)) / pow(Zp1,4); // Formula as in the paper. But I implement it using the ratio of IonizPotentials. See Golovanov 2020 paper
    // The rate computed below is in  AtomicUnits. E-field received as argument is in AtomicUnits (see IonizationTunnelBSI::operator()(arg1, arg2, ...) function)
    double Z = Zp1 - 1;
    double IH = 13.598434005136; // IP of atomic Hydrogen H, in eV.
    double Ii = IonizationTables::ionization_energy(atomic_number_, Z);
    double ratio_of_IPs = IH / Ii;
    
    double BSI_rate_quadratic = 2.4 * (pow(E,2)) * pow(ratio_of_IPs,2); // E is in atomic units here, rate is calculated in atomic units
    // return BSI_rate_quadratic; // Returns the rate in ATOMIC UNITS (i.e. normalized to atomic_unit_of_frequency)
    return (BSI_rate_quadratic * au_to_w0); // in SMILEI UNITS
}

// Derivatives of BSI Quadratic rate:
// double BSI_rate_quadratic_function_der_wrtE(unsigned int Zp1, double E, unsigned int atomic_number_) { // simple derivative of BSI quadratic rate wrt E-field 
//     // double BSI_rate_quadratic_der_wrtE = 2.4 * 2 * E / pow(Zp1,4); // The derivate of the formula as it's appearing in the paper by Bauer-Mulser BM.
//     double Z = Zp1 - 1;
//     double IH = 13.598434005136;
//     double Ii = IonizationTables::ionization_energy(atomic_number_, Z);
//     double ratio_of_IPs = IH / Ii;

//     double BSI_rate_quadratic_der_wrtE = 2.4 * 2 * E *  pow(ratio_of_IPs,2);
//     return (BSI_rate_quadratic_der_wrtE * au_to_w0); // Returned in SMILEI UNITS.
//     // return BSI_rate_quadratic_der_wrtE; // Returned in Atomic Units.
// }

// double BSI_rate_quadratic_function_loglogder_wrtE(unsigned int Zp1, double E, unsigned int atomic_number_) {
//     double BSI_rate_quadratic = BSI_rate_quadratic_function(Zp1,E, atomic_number_);
//     double BSI_rate_quadratic_der_wrtE = BSI_rate_quadratic_function_der_wrtE(Zp1,E, atomic_number_);

//     double BSI_rate_quadratic_loglogder_wrtE = (E/BSI_rate_quadratic) * BSI_rate_quadratic_der_wrtE; // log-log derivative is: (x/y) * (dy/dx) where in our case: y = Rate, x = E-field.
//     return BSI_rate_quadratic_loglogder_wrtE;
// }



// 2) BSI-LINEAR dependence on local E-field <---> Kostyukov+Golovanov: Field ionization in short and extremely intense laser pulses, 2018
double BSI_rate_linear_function(unsigned int Zp1, double E, unsigned int atomic_number_, double au_to_w0) { // E-field E from arguments has to be in AtomicUnits. 
//     double BSI_rate_linear = 0.8 * (E / Zp1) ; // 1/Zp1 = sqrt(IH[eV] / Ii[ev])
        double Z = Zp1 - 1;
        double IH = 13.598434005136; // IP of atomic Hydrogen H, in eV.
        double Ii = IonizationTables::ionization_energy(atomic_number_, Z); // Ii shall be the IP of atom/ion. Instr. on RSH of equal sign returns IP of element with atomic no = atomic_number_ and charge state Z (Z=0 is neutral atom).
        double ratio_of_IPs = IH / Ii;
        
        double BSI_rate_linear = 0.8 * E * pow(ratio_of_IPs, 0.5); // E is in atomic units, rate is computed in Atomic units.
        // return BSI_rate_linear; // Returns BSI linear rate in AtomicUnits (i.e. normalized to atomic_unit_of_frequency)
        return (BSI_rate_linear * au_to_w0); // in SMILEI units.
}

// Derivatives of BSI Linear rate:
// double BSI_rate_linear_function_der_wrtE(unsigned int Zp1, unsigned int atomic_number_) { // E-field E from arguments has to be in AtomicUnits
//     double Z = Zp1 - 1;
//     double IH = 13.598434005136; // IP of atomic Hydrogen H, in eV.
//     double Ii = IonizationTables::ionization_energy(atomic_number_, Z); // Ii shall be the IP of atom/ion. Instr. on RSH of equal sign returns IP of element with atomic no = atomic_number_ and charge state Z (Z=0 is neutral atom).
//     double ratio_of_IPs = IH / Ii;
    
//     double BSI_rate_linear_der_wrtE = 0.8 * pow(ratio_of_IPs, 0.5);
//     return BSI_rate_linear_der_wrtE * au_to_w0;
//     // return BSI_rate_linear_der_wrtE;
// } 

// double BSI_rate_linear_function_loglogder_wrtE(unsigned int Zp1, double E, unsigned int atomic_number_) {
//     double BSI_linear_rate = BSI_rate_linear_function(Zp1, E, atomic_number_);
//     double BSI_linear_rate_der_wrtE = BSI_rate_linear_function_der_wrtE(Zp1, atomic_number_);

//     double BSI_linear_rate_loglogder_wrtE = (E / BSI_linear_rate) * BSI_linear_rate_der_wrtE;

//     return BSI_linear_rate_loglogder_wrtE;
// }
