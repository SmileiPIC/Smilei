#ifndef PARTICLE_H
#define PARTICLE_H

#include <cstring>
#include <iostream>
#include <vector>

class Particles;

class Particle
{
public:
    Particle() {};
    Particle( Particles &parts, int iPart );
    
    ~Particle() {};
    
    friend std::ostream &operator<<( std::ostream &os, const Particle &part );
    
private:
    //! array containing the particle position
    std::vector<double> Position;
    //! array containing the particle former (old) positions
    std::vector<double> Position_old;
    //! array containing the particle moments
    std::vector<double>  Momentum;
    //! containing the particle weight: equivalent to a charge density
    double Weight;
    //! containing the particle quantum parameter
    double Chi;
    //! containing the particle optical depth
    double Tau;
    //! charge state of the particle (multiples of e>0)
    short Charge;
    //! Id of the particle
    uint64_t Id;
};

#endif
