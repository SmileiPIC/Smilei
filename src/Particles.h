#ifndef PARTICLES_H
#define PARTICLES_H

#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include "Tools.h"

class PicParams;

//----------------------------------------------------------------------------------------------------------------------
//! Particle class: holds the basic properties of a particle
//----------------------------------------------------------------------------------------------------------------------
class Particles {
public:

    //! Constructor for Particle
    Particles();

    //! Destructor for Particle
    virtual ~Particles();

    //! Create nParticles null particles of nDim size
    void initialize( int nParticles, int nDim );
    //! Set capacity of Particles vectors
    void reserve( unsigned int n_part_max, int nDim );

    // Reset Particles vectors
    void clear();

    // Get number of particules
    inline int size() const {
        return Weight.size();
    }
    // Get dimension of particules
    inline int dimension() const {
        return Position.size();
    }

    // Copy particle iPart at the end of dest_parts
    void cp_particle(int iPart, Particles &dest_parts );
    // Copy particle iPart at dest_id in dest_parts
    void cp_particle(int ipart, Particles &dest_parts, int dest_id );

    // Suppress particle iPart
    void erase_particle(int iPart );

    // Print parameters of particle iPart
    void print(int iPart);

    // Exchange particles part1 & part2 memory location
    void swap_part(int part1,int part2);

    // Move iPart at the end of vectors
    void push_to_end(int iPart );

    // Create new particle
    void create_particle();
    // Create nParticles new particles
    void create_particles(int nParticles);

    //! Method used to get the Particle position
    inline double  position( int idim, int ipart ) const {
        return Position[idim][ipart];
    }
    //! Method used to set a new value to the Particle former position
    inline double& position( int idim, int ipart )       {
        return Position[idim][ipart];
    }

    //! Method used to get the Particle position
    inline double  position_old( int idim, int ipart ) const {
        return Position_old[idim][ipart];
    }
    //! Method used to set a new value to the Particle former position
    inline double& position_old( int idim, int ipart )       {
        return Position_old[idim][ipart];
    }


    //! Method used to get the Particle momentum
    inline double  momentum( int idim, int ipart ) const {
        return Momentum[idim][ipart];
    }
    //! Method used to set a new value to the Particle momentum
    inline double& momentum( int idim, int ipart )       {
        return Momentum[idim][ipart];
    }

    //! Method used to get the Particle weight
    inline double  weight(int ipart) const {
        return Weight[ipart];
    }
    //! Method used to set a new value to the Particle weight
    inline double& weight(int ipart)       {
        return Weight[ipart];
    }

    //! Method used to get the Particle charge
    inline short  charge(int ipart) const {
        return Charge[ipart];
    }
    //! Method used to set a new value to the Particle weight
    inline short& charge(int ipart)       {
        return Charge[ipart];
    }

    //! Method used to get the Particle Lorentz factor
    inline double lor_fac(int ipart) {
        return sqrt(1+pow(momentum(0,ipart),2)+pow(momentum(1,ipart),2)+pow(momentum(2,ipart),2));
    }

private:
    //! array containing the particle position
    std::vector< std::vector<double> > Position;

    //! array containing the particle former (old) positions
    std::vector< std::vector<double> >Position_old;

    //! array containing the particle moments
    std::vector< std::vector<double> >  Momentum;

    //! containing the particle weight: equivalent to a charge density
    std::vector<double> Weight;

    //! charge state of the particle (multiples of e>0)
    std::vector<short> Charge;

};

#endif
