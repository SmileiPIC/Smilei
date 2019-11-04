#ifndef IONIZATIONTUNNELENVELOPEAVERAGED_H
#define IONIZATIONTUNNELENVELOPEAVERAGED_H

#include <cmath>

#include <vector>

#include "Ionization.h"
#include "Tools.h"
//#include <boost/math/special_functions/bessel.hpp>

class Particles;

//! calculate the particle tunnel ionization
class IonizationTunnelEnvelopeAveraged : public Ionization
{

public:
    //! Constructor for IonizationTunnelEnvelope: with no input argument
    IonizationTunnelEnvelopeAveraged( Params &params, Species *species );
    
    //! apply the Tunnel Ionization model to the species (with ionization current)
    void operator()( Particles *, unsigned int, unsigned int, std::vector<double> *, Patch *, Projector *, int ipart_ref = 0 ) override;
    //! method for envelope ionization
    void envelopeIonization( Particles *, unsigned int, unsigned int, std::vector<double> *, std::vector<double> *, std::vector<double> *, Patch *, Projector *, int ipart_ref = 0 ) override;

    inline double correction_factor_general_polarization(double x){
        return exp(-x);//*boost::math::cyl_bessel_i(0, x); 
    };

    double ellipticity,cos_phi,sin_phi,ellipticity_factor1,ellipticity_factor2;

private:
    unsigned int atomic_number_;
    std::vector<double> Potential;
    std::vector<double> Azimuthal_quantum_number;
    
    double one_third;
    std::vector<double> alpha_tunnel, beta_tunnel, gamma_tunnel; //,Ip_times2_to_minus3ov4;
};


#endif
