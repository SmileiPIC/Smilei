#ifndef PARTICLES_H
#define PARTICLES_H

#include <cmath>

#include <iostream>
#include <fstream>
#include <vector>

#include "Tools.h"
#include "TimeSelection.h"

class Particle;

class Params;
class Patch;



//----------------------------------------------------------------------------------------------------------------------
//! Particle class: holds the basic properties of a particle
//----------------------------------------------------------------------------------------------------------------------
class Particles
{
public:

    //! Constructor for Particle
    Particles();

    //! Destructor for Particle
    ~Particles() {}

    //! Create nParticles null particles of nDim size
    void initialize( unsigned int nParticles, unsigned int nDim );

    //! Create nParticles null particles of nDim size
    void initialize( unsigned int nParticles, Particles &part );

    //! Set capacity of Particles vectors
    void reserve( unsigned int n_part_max, unsigned int nDim );
    
    //! Initialize like another particle, but only reserve space
    void initialize_reserve( unsigned int n_part_max, Particles &part );

    //! Resize Particles vectors
    void resize( unsigned int nParticles, unsigned int nDim );

    //! Resize Particles vectors
    void resize( unsigned int nParticles);

    //! Remove extra capacity of Particles vectors
    void shrink_to_fit();

    //! Reset Particles vectors
    void clear();

    //! Get number of particules
    inline unsigned int size() const
    {
        return Weight.size();
    }

    //! Get number of particules
    inline unsigned int capacity() const
    {
        return Weight.capacity();
    }

    //! Get dimension of particules
    inline unsigned int dimension() const
    {
        return Position.size();
    }

    //! Copy particle iPart at the end of dest_parts
    void cp_particle( unsigned int iPart, Particles &dest_parts );

    void cp_particle( unsigned int iPart );


    //! Insert nPart particles starting at ipart to dest_id in dest_parts
    void cp_particles( unsigned int iPart, unsigned int nPart, Particles &dest_parts, int dest_id );
    //! Insert particle iPart at dest_id in dest_parts
    void cp_particle( unsigned int ipart, Particles &dest_parts, int dest_id );
    
    //! Copy particle iPart at the end of dest_parts -- safe
    void cp_particle_safe( unsigned int ipart, Particles &dest_parts );
    
    //! Suppress particle iPart
    void erase_particle( unsigned int iPart );
    //! Suppress nPart particles from iPart
    void erase_particle( unsigned int iPart, unsigned int nPart );

    //! Suppress all particles from iPart to the end of particle array
    void erase_particle_trail( unsigned int iPart );

    //! Print parameters of particle iPart
    void print( unsigned int iPart );

    friend std::ostream &operator << ( std::ostream &, const Particles &particle );

    //! Exchange particles part1 & part2 memory location
    void swap_part( unsigned int part1, unsigned int part2 );
    void swap_parts( std::vector<unsigned int> parts );
    void translate_parts( std::vector<unsigned int> parts );
    void swap_part3( unsigned int part1, unsigned int part2, unsigned int part3 );
    void swap_part4( unsigned int part1, unsigned int part2, unsigned int part3, unsigned int part4 );

    //! Exchange particles part1 & part2 memory location
    void swap_part( unsigned int part1, unsigned int part2, unsigned int N );

    //! Overwrite particle part1 into part2 memory location. Erasing part2
    void overwrite_part( unsigned int part1, unsigned int part2 );

    //! Overwrite particle part1->part1+N into part2->part2+N memory location. Erasing part2->part2+N
    void overwrite_part( unsigned int part1, unsigned int part2, unsigned int N );

    //! Overwrite particle part1->part1+N into part2->part2+N of dest_parts memory location. Erasing part2->part2+N
    void overwrite_part( unsigned int part1, Particles &dest_parts, unsigned int part2, unsigned int N );

    //! Overwrite particle part1 into part2 of dest_parts memory location. Erasing part2
    void overwrite_part( unsigned int part1, Particles &dest_parts, unsigned int part2 );


    //! Move iPart at the end of vectors
    void push_to_end( unsigned int iPart );

    //! Create new particle
    void create_particle();

    //! Create nParticles new particles
    void create_particles( int nAdditionalParticles );

    //! Compress the particles vectors according to the provided mask
    //! between istart and iend
    void compressParticles( int istart, int iend, std::vector <int> & mask );

    //! Compress the particles vectors using cell_keys as a mask
    //! between istart and iend
    void compressParticles( int istart, int iend);

    //! Test if ipart is in the local patch
    bool is_part_in_domain( unsigned int ipart, Patch *patch );

    //! Method used to get the Particle position
    inline double  position( unsigned int idim, unsigned int ipart ) const
    {
        return Position[idim][ipart];
    }
    //! Method used to set a new value to the Particle former position
    inline double &position( unsigned int idim, unsigned int ipart )
    {
        return Position[idim][ipart];
    }

    //! Method used to get the Particle position
    inline double distance2_to_axis( unsigned int ipart ) const
    {
        return Position[1][ipart] * Position[1][ipart] + Position[2][ipart] * Position[2][ipart];
    }

    //! Method used to get the Particle position
    inline double  position_old( unsigned int idim, unsigned int ipart ) const
    {
        return Position_old[idim][ipart];
    }
    //! Method used to set a new value to the Particle former position
    inline double &position_old( unsigned int idim, unsigned int ipart )
    {
        return Position_old[idim][ipart];
    }

    //! Method used to get the list of Particle position
    inline std::vector<double>  position( unsigned int idim ) const
    {
        return Position[idim];
    }

    //! Method used to get the Particle momentum
    inline double  momentum( unsigned int idim, unsigned int ipart ) const
    {
        return Momentum[idim][ipart];
    }
    //! Method used to set a new value to the Particle momentum
    inline double &momentum( unsigned int idim, unsigned int ipart )
    {
        return Momentum[idim][ipart];
    }
    //! Method used to get the Particle momentum
    inline std::vector<double>  momentum( unsigned int idim ) const
    {
        return Momentum[idim];
    }

    //! Method used to get the Particle weight
    inline double  weight( unsigned int ipart ) const
    {
        return Weight[ipart];
    }
    //! Method used to set a new value to the Particle weight
    inline double &weight( unsigned int ipart )
    {
        return Weight[ipart];
    }
    //! Method used to get the Particle weight
    inline std::vector<double>  weight() const
    {
        return Weight;
    }

    //! Method used to get the Particle charge
    inline short  charge( unsigned int ipart ) const
    {
        return Charge[ipart];
    }
    //! Method used to set a new value to the Particle charge
    inline short &charge( unsigned int ipart )
    {
        return Charge[ipart];
    }
    //! Method used to get the list of Particle charges
    inline std::vector<short>  charge() const
    {
        return Charge;
    }


    //! Method used to get the Particle Lorentz factor
    inline double lor_fac( unsigned int ipart )
    {
        return sqrt( 1.+pow( momentum( 0, ipart ), 2 )+pow( momentum( 1, ipart ), 2 )+pow( momentum( 2, ipart ), 2 ) );
    }

    //! Method used to get the inverse Particle Lorentz factor
    inline double inv_lor_fac( unsigned int ipart )
    {
        return 1./sqrt( 1.+pow( momentum( 0, ipart ), 2 )+pow( momentum( 1, ipart ), 2 )+pow( momentum( 2, ipart ), 2 ) );
    }

    //! Method used to get the momentum norm which is also the normalized photon energy
    inline double momentum_norm( unsigned int ipart )
    {
        return sqrt( pow( momentum( 0, ipart ), 2 )+pow( momentum( 1, ipart ), 2 )+pow( momentum( 2, ipart ), 2 ) );
    }

    //! Partiles properties, respect type order : all double, all short, all unsigned int

    //! array containing the particle position
    std::vector< std::vector<double> > Position;

    //! array containing the particle former (old) positions
    std::vector< std::vector<double> >Position_old;

    //! array containing the particle moments
    std::vector< std::vector<double> >  Momentum;

    //! containing the particle weight: equivalent to a charge density
    std::vector<double> Weight;

    //! containing the particle quantum parameter
    std::vector<double> Chi;

    //! charge state of the particle (multiples of e>0)
    std::vector<short> Charge;

    //! Id of the particle
    std::vector<uint64_t> Id;

    // Discontinuous radiation losses

    //! Incremental optical depth for
    //! the Monte-Carlo process
    std::vector<double> Tau;

    //! cell_keys of the particle
    std::vector<int> cell_keys;

    // TEST PARTICLE PARAMETERS
    bool is_test;

    //! True if tracking the particles
    bool tracked;

    void resetIds()
    {
        unsigned int s = Id.size();
        for( unsigned int iPart=0; iPart<s; iPart++ ) {
            Id[iPart] = 0;
        }
    }

    //! Method used to get the Particle Id
    inline uint64_t id( unsigned int ipart ) const
    {
        DEBUG( ipart << " of " << Id.size() );
        return Id[ipart];
    }
    //! Method used to set the Particle Id
    inline uint64_t &id( unsigned int ipart )
    {
        return Id[ipart];
    }
    //! Method used to get the Particle Ids
    inline std::vector<uint64_t> id() const
    {
        return Id;
    }
    void sortById();

    //! Quantum parameter for particles that are submitted
    //! to a radiation reaction force (CED or QED)
    bool isQuantumParameter;

    //! Parameters for particles that are submitted to a
    //! Monte-Carlo process such as:
    //! - discontinuous radiation reaction force
    bool isMonteCarlo;

    //! Method used to get the Particle chi factor
    inline double  chi( unsigned int ipart ) const
    {
        return Chi[ipart];
    }
    //! Method used to set a new value to the Particle chi factor
    inline double &chi( unsigned int ipart )
    {
        return Chi[ipart];
    }
    //! Method used to get the Particle chi factor
    inline std::vector<double>  chi() const
    {
        return Chi;
    }

    //! Method used to get the Particle optical depth
    inline double  tau( unsigned int ipart ) const
    {
        return Tau[ipart];
    }
    //! Method used to set a new value to
    //! the Particle optical depth
    inline double &tau( unsigned int ipart )
    {
        return Tau[ipart];
    }
    //! Method used to get the Particle optical depth
    inline std::vector<double>  tau() const
    {
        return Tau;
    }


    std::vector< std::vector<double  >*> double_prop;
    std::vector< std::vector<short   >*> short_prop;
    std::vector< std::vector<uint64_t>*> uint64_prop;


#ifdef __DEBUG
    bool test_move( int iPartStart, int iPartEnd, Params &params );

    inline double dist2( unsigned int iPart )
    {
        double dist( 0. );
        for( unsigned int iDim = 0 ; iDim < Position.size() ; iDim++ ) {
            double delta = position( iDim, iPart )-position_old( iDim, iPart );
            dist += delta*delta;
        }
        return dist;
    }
    inline double dist( unsigned int iPart, unsigned int iDim )
    {
        double delta = std::abs( position( iDim, iPart )-position_old( iDim, iPart ) );
        return delta;
    }
#endif

    Particle operator()( unsigned int iPart );

    //! Methods to obtain any property, given its index in the arrays double_prop, uint64_prop, or short_prop
    void getProperty( unsigned int iprop, std::vector<uint64_t> *&prop )
    {
        prop = uint64_prop[iprop];
    }
    void getProperty( unsigned int iprop, std::vector<short> *&prop )
    {
        prop = short_prop[iprop];
    }
    void getProperty( unsigned int iprop, std::vector<double> *&prop )
    {
        prop = double_prop[iprop];
    }

private:

};



#endif
