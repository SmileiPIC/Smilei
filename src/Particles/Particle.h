#ifndef PARTICLE_H
#define PARTICLE_H

#include <cstring>
#include <iostream>
#include <vector>
#include <cstdint>

class Particles;

class Particle
{
public:
    Particle() {};
    Particle( Particles &parts, int iPart );
    
    ~Particle() {};
    
    friend std::ostream &operator<<( std::ostream &os, const Particle &part );
    
private:
    //! particle position
    std::vector<double> Position;
    //! particle former (old) position
    std::vector<double> Position_old;
    //! particle momentum
    std::vector<double>  Momentum;
    //! particle weight: equivalent to a density normalized to the number of macro-particles per cell
    double Weight;
    //! particle quantum parameter
    double Chi;
    //! particle optical depth for Monte-Carlo processes
    double Tau;
    //! particle charge state
    short Charge;
    //! particle Id
    uint64_t Id;
};

#endif
