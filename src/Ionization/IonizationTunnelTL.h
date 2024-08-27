#ifndef IONIZATIONTUNNELTL_H
#define IONIZATIONTUNNELTL_H

#include <cmath>

#include <vector>

#include "Ionization.h"
#include "Tools.h"

class Particles;

//! calculate the particle tunnel ionization with Tong-Lin model
class IonizationTunnelTL : public Ionization
{

public:
    //! Constructor for IonizationTunnelTL: with no input argument
    IonizationTunnelTL( Params &params, Species *species );
    
    //! apply the Tunnel Ionization model to the species (with ionization current)
    void operator()( Particles *, unsigned int, unsigned int, std::vector<double> *, Patch *, Projector *, int ipart_ref = 0 ) override;
    
private:
    unsigned int atomic_number_;
    std::vector<double> Potential;
    std::vector<double> Azimuthal_quantum_number;
    
    double one_third;
    double ionization_tl_parameter_;
    std::vector<double> alpha_tunnel, beta_tunnel, gamma_tunnel, lambda_tunnel;
};


#endif
