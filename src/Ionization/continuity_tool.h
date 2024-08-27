#ifndef CONTINUITYTOOL_H
#define CONTINUITYTOOL_H

// #include "Ionization.h"

double BSI_rate_quadratic_function(unsigned int Zp1, double E, unsigned int atomic_number_, double au_to_w0);
//double BSI_rate_quadratic_function_der_wrtE(unsigned int Zp1, double E, unsigned int atomic_number_);
//double BSI_rate_quadratic_function_loglogder_wrtE(unsigned int Zp1, double E, unsigned int atomic_number_);

double BSI_rate_linear_function(unsigned int Zp1, double E, unsigned int atomic_number_, double au_to_w0);
//double BSI_rate_linear_function_der_wrtE(unsigned int Zp1, unsigned int atomic_number_);
//double BSI_rate_linear_function_loglogder_wrtE(unsigned int Zp1, double E, unsigned int atomic_number_);

double Tunnel_rate_function(unsigned int Z, double E, double alpha, double beta, double gamma);
//double Tunnel_rate_function_der_wrtE(unsigned int Z, double E, double alpha, double beta, double gamma);
//double Tunnel_rate_function_loglogder_wrtE(unsigned int Z, double E, double alpha, double beta, double gamma);


int continuity_tool(unsigned int Zp1, double E, double alpha, double beta, double gamma, double E_cr, double Potential, unsigned int atomic_number_, double au_to_w0);

extern double au_to_eV; // au_to_eV will be defined in continuity_tool.cpp (and defined only there, only once)
                        // will need to #include "continuity_tool.h" in each .cpp that uses au_to_eV.
#endif