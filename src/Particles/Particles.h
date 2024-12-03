// -----------------------------------------------------------------------------
//
//! \file Particles.cpp
//
//! \brief contains the Particles class description
//
//! The Particles Class is the main data structure for handling particle list.
//! It contains the main particles properties:
//! - positions
//! - momentums
//! - charge
//! - weight
//! - quantum parameter (chi) for QED effects
//! - optical depth for Monte-Carlo processes
//! - tag id for tracked particles
//
//! The class also contains many functions to manage particles.
// -----------------------------------------------------------------------------

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

struct InterpolatedFields {
    //! Tells the way each interpolated field is treated: 0 = not kept, 1 = kept, 2 = accumulated
    std::vector<int> mode_;
    //! arrays of fields interpolated on the particle positions. The order is Ex, Ey, Ez, Bx, By, Bz, Wx, Wy, Wz
    std::vector<std::vector<double>> F_;
};


//----------------------------------------------------------------------------------------------------------------------
//! Particle class: holds the basic properties of a particle
//----------------------------------------------------------------------------------------------------------------------
class Particles
{
public:

    //! Constructor for Particle
    Particles();

    //! Destructor for Particle
    virtual ~Particles();

    //! Create nParticles null particles of nDim size
    void initialize( unsigned int nParticles, unsigned int nDim, bool keep_position_old );

    //! Create nParticles null particles of nDim size
    void initialize( unsigned int nParticles, Particles &part );

    //! Set capacity of Particles vectors and change dimensionality
    void reserve( unsigned int reserved_particles, unsigned int nDim, bool keep_position_old = false);

    //! Set capacity of Particles vectors and keep dimensionality
    void reserve( unsigned int reserved_particles );

    //! Initialize like Particles object part with 0 particles and reserve space for n_part_max particles
    void initializeReserve( unsigned int n_part_max, Particles &part );

    //! Resize Particle vectors and change dimensionality according to nDim
    void resize( unsigned int nParticles, unsigned int nDim, bool keep_position_old );

    //! Resize Particles vectors
    void resize( unsigned int nParticles);

    //! Resize the cell_keys vector
    void resizeCellKeys(unsigned int nParticles);

    //! Remove extra capacity of Particles vectors
    //! params [in] compute_cell_keys: if true, cell_keys is affected (default is false)
    void shrinkToFit(const bool compute_cell_keys = false);

    //! Reset Particles vectors
    //! params [in] compute_cell_keys: if true, cell_keys is affected (default is false)
    void clear(const bool compute_cell_keys = false);

    //! Get number of particles
    inline unsigned int numberOfParticles() const
    {
        // If the notion of bin is not used, the vector size is the number of Particles
        if (last_index.size() == 0) {
            //ERROR("Particles object tried to use `numberOfParticles` but `last_index` is not initialized.")
            return Weight.size();
        }
        return last_index.back();
    }

    //! Get vector size on CPU
    inline unsigned int size() const
    {
        return Weight.size();
    }

    //! Get vector size on CPU
    inline unsigned int hostVectorSize() const
    {
        return Weight.size();
    }

    //! Get number of particles
    inline unsigned int capacity() const
    {
        return Weight.capacity();
    }

    //! Get dimension of particles
    inline unsigned int dimension() const
    {
        return Position.size();
    }

    //! Get dimension of particles
    inline unsigned int numberOfBins() const
    {
        return first_index.size();
    }

    //! Tells if old positions are kept (true) or not
    inline bool keepOldPositions() const
    {
        return Position_old.size() > 0;
    }

    //! Copy particle iPart at the end of dest_parts
    void copyParticle( unsigned int iPart, Particles &dest_parts );
    //! Copy particle iPart at the end of the current array
    void copyParticle( unsigned int iPart );
    //! Insert particle iPart at dest_id in dest_parts
    void copyParticle( unsigned int ipart, Particles &dest_parts, int dest_id );

    //! Insert nPart particles starting at ipart to dest_id in dest_parts
    void copyParticles( unsigned int iPart, unsigned int nPart, Particles &dest_parts, int dest_id );
    //! Transfer particles indexed by array indices to dest_id in dest_parts
    void copyParticles( std::vector<size_t> indices, Particles &dest_parts, int dest_id );

    //! Make a new particle at the position of another
    void makeParticleAt( Particles &source_particles, unsigned int ipart, double w, short q=0., double px=0., double py=0., double pz=0. );

    //! Suppress particle iPart
    void eraseParticle( unsigned int iPart, bool compute_cell_keys = false );
    //! Suppress nPart particles from iPart
    void eraseParticle( unsigned int iPart, unsigned int nPart, bool compute_cell_keys = false );
    //! Suppress indexed particles
    void eraseParticles( std::vector<size_t> indices );

    //! Suppress all particles from iPart to the end of particle array
    void eraseParticleTrail( unsigned int iPart, bool compute_cell_keys = false );

    //! Print parameters of particle iPart
    void print( unsigned int iPart );

    friend std::ostream &operator << ( std::ostream &, const Particles &particle );

    //! Exchange particles part1 & part2 memory location
    void swapParticle( unsigned int part1, unsigned int part2 );
    void swapParticles( std::vector<unsigned int> parts );
    void translateParticles( std::vector<unsigned int> parts );
    void swapParticle3( unsigned int part1, unsigned int part2, unsigned int part3 );
    void swapParticle4( unsigned int part1, unsigned int part2, unsigned int part3, unsigned int part4 );

    //! Exchange particles part1 & part2 memory location
    void swapParticle( unsigned int part1, unsigned int part2, unsigned int N );

    //! Overwrite particle part1 into part2 memory location. Erasing part2
    //! Warning: do not update first_index and last_index
    void overwriteParticle( unsigned int part1, unsigned int part2, bool compute_cell_keys = false  );

    //! Overwrite particle part1->part1+N into part2->part2+N memory location. Erasing part2->part2+N
    //! Warning: do not update first_index and last_index
    void overwriteParticle( unsigned int part1, unsigned int part2, unsigned int N, bool compute_cell_keys = false );

    //! Overwrite particle part1->part1+N into part2->part2+N of dest_parts memory location. Erasing part2->part2+N
    //! Warning: do not update first_index and last_index
    void overwriteParticle( unsigned int part1, Particles &dest_parts, unsigned int part2, unsigned int N );

    //! Overwrite particle part1 into part2 of dest_parts memory location. Erasing part2
    //! Warning: do not update first_index and last_index
    void overwriteParticle( unsigned int part1, Particles &dest_parts, unsigned int part2 );

    //! Create new particle
    void createParticle();

    //! Create n_additional_particles new particles
    virtual void createParticles( int n_additional_particles );

    //! Create n_additional_particles new particles at position pstart in the particles data structure
    void createParticles( int n_additional_particles, int pstart );

    //! Move ipart at new_pos in the particles data structure
    void moveParticles( int iPart, int new_pos );

    //! Remove and compress the particles vectors according to the provided mask
    //! between istart and iend
    void eraseParticlesWithMask( int istart, int iend, std::vector <int> & mask );

    //! Remove and compress the particles vectors using cell_keys as a mask
    //! between istart and iend
    void eraseParticlesWithMask( int istart, int iend);

    //! This method erases particles according to the provided mask
    //! between istart and iend
    // void eraseParticlesWithMask( int istart, int iend, vector <bool> & to_be_erased);

    //! This method eliminates the space between the bins
    //! (presence of empty particles beteen the bins)
    void compress(bool compute_cell_keys = false);

    //! Sum the vectors
    void sum(int ibin_min, int ibin_max);

    //! Test if ipart is in the local patch
    bool isParticleInDomain( unsigned int ipart, Patch *patch );

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
    inline double distance2ToAxis( unsigned int ipart ) const
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
    inline  double LorentzFactor( unsigned int ipart )
    {
        return sqrt( 1. + momentum( 0, ipart ) * momentum( 0, ipart ) + momentum( 1, ipart ) * momentum( 1, ipart ) + momentum( 2, ipart ) * momentum( 2, ipart )  );
    }

    //! Method used to get the inverse Particle Lorentz factor
    inline double inverseLorentzFactor( unsigned int ipart )
    {
        return 1./sqrt( 1.+ momentum( 0, ipart )*momentum( 0, ipart ) + momentum( 1, ipart )*momentum( 1, ipart ) + momentum( 2, ipart )*momentum( 2, ipart ) );
    }

    //! Method used to get the momentum norm which is also the normalized photon energy
    inline double momentumNorm( unsigned int ipart )
    {
        return sqrt( momentum( 0, ipart ) * momentum( 0, ipart ) + momentum( 1, ipart ) * momentum( 1, ipart ) + momentum( 2, ipart ) * momentum( 2, ipart ) );
    }

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

    //! Method to keep the positions for the next timesteps
    void savePositions();

    std::vector< std::vector<double  >*> double_prop_;
    std::vector< std::vector<short   >*> short_prop_;
    std::vector< std::vector<uint64_t>*> uint64_prop_;

#ifdef __DEBUG
    bool testMove( int iPartStart, int iPartEnd, Params &params );

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
    
    void prepareInterpolatedFields( std::vector<std::vector<double>> &pold, size_t start, size_t n );
    void copyInterpolatedFields( double *Ebuffer, double *Bbuffer, std::vector<std::vector<double>> &pold, size_t start, size_t n, size_t buffer_size, double mass_ );

    //! Methods to obtain any property, given its index in the arrays double_prop_, uint64_prop_, or short_prop_
    void getProperty( size_t iprop, std::vector<uint64_t> *&prop )
    {
        prop = uint64_prop_[iprop];
    }
    void getProperty( size_t iprop, std::vector<short> *&prop )
    {
        prop = short_prop_[iprop];
    }
    void getProperty( size_t iprop, std::vector<double> *&prop )
    {
        prop = double_prop_[iprop];
    }

    //! Make a copy of the host particles and does some computation like
    //! bin discovery.
    //!
    virtual void initializeDataOnDevice();
    virtual void initializeIDsOnDevice();
    virtual void copyFromHostToDevice();
    virtual void copyFromDeviceToHost( bool copy_keys = false );

    //! Return the pointer toward the Position[idim] vector
    virtual double* getPtrPosition( int idim ) {
        return ((std::size_t)idim < Position.size()) ? Position[idim].data() : nullptr;
    };
    //! Return the pointer toward the Position_old[idim] vector
    virtual double* getPtrPositionOld( int idim ) {
        return ((std::size_t)idim < Position_old.size()) ? Position_old[idim].data() : nullptr;
    };
    //! Return the pointer toward the Momentum[idim] vector
    virtual double* getPtrMomentum( int idim ) {
        return ((std::size_t)idim < Momentum.size()) ? Momentum[idim].data() : nullptr;
    };
    virtual double* getPtrWeight() {
        return &(Weight[0]);
    };
    virtual double* getPtrChi() {
        return (has_quantum_parameter ? Chi.data() : nullptr);
    };
    virtual short* getPtrCharge() {
        return &(Charge[0]);
    };
    virtual uint64_t* getPtrId() {
        return &(Id[0]);
    };
    virtual double* getPtrTau() {
        return (has_Monte_Carlo_process ? Tau.data() : nullptr);
    };
    virtual int* getPtrCellKeys() {
        return &(cell_keys[0]);
    };


    // --------------------------------------------------------------------------------------------
    // Accelerator specific virtual functions

    // -----------------------------------------------------------------------------
    //! Extract particles leaving the box to buffers
    // -----------------------------------------------------------------------------
    void copyLeavingParticlesToBuffers( const std::vector<bool> copy, const std::vector<Particles*> buffer );
    virtual void copyLeavingParticlesToBuffer( Particles* buffer );

    // -----------------------------------------------------------------------------
    //! Erase particles leaving the patch object on device
    // -----------------------------------------------------------------------------
    virtual int eraseLeavingParticles();

    // -----------------------------------------------------------------------------
    //! Resize & Copy particles from particles_to_inject to the end of the vectors
    virtual int addParticles( Particles* particles_to_inject  );
    
    //! Implementation of a somewhat efficient particle injection, sorting
    //! (including removing leaving particles) and binning for GPU if
    //! available for the configuration of offloading technology
    //! (OpenMP/OpenACC) and implemented for the space dimension of the
    //! simulation. Else, fallback to the old naive, plain memcpy
    //! implementation.
    //!
    //! last_index is modified appropriately.
    //!
    virtual void importAndSortParticles( Particles *particles_to_inject );

    //! Returns the capacity of the vectors representing the particle
    //! components. One can resize up to deviceCapacity without triggering
    //! a reallocation and potentially invaliding the iterators.
    virtual unsigned int deviceCapacity() const;

    //! Get number of particles on device
    virtual unsigned int deviceSize() const {
        ERROR( "deviceSize is a feature only available for accelerator device" );
        return 0;
    }

    virtual void setHostBinIndex();

    // ---------------------------------------------------------------------------------------
    // Parameters
    // partiles properties, respect type order : all double, all short, all unsigned int

    //! array of particle positions
    std::vector< std::vector<double> > Position;

    //! array of particle former (old) positions
    std::vector< std::vector<double> >Position_old;

    //! array of particle momenta
    std::vector< std::vector<double> >  Momentum;

    //! array of particle weights: equivalent to a density normalized to the number of macro-particles per cell
    std::vector<double> Weight;

    //! array of particle quantum parameters
    std::vector<double> Chi;

    //! array of optical depths for the Monte-Carlo process
    std::vector<double> Tau;

    //! array of particle charges
    std::vector<short> Charge;

    //! array of particle IDs
    std::vector<uint64_t> Id;
    
    //! arrays of fields interpolated at particle positions
    InterpolatedFields * interpolated_fields_;
    
    //! array of particle cell keys (for sorting per cell)
    std::vector<int> cell_keys;

    // TEST PARTICLE PARAMETERS
    bool is_test;

    //! True if tracking the particles
    bool tracked;

    //! Indices of the first particles of each bin (or cells) in the Particles object
    std::vector<int> first_index;

    //! Indexes of the last particles + 1 in each bin (or cells) in the Particles object
    std::vector<int> last_index;

    //! Quantum parameter for particles that are submitted
    //! to a radiation reaction force (CED or QED)
    bool has_quantum_parameter;

    //! Parameters for particles that are submitted to a
    //! Monte-Carlo process such as:
    //! - discontinuous radiation reaction force
    bool has_Monte_Carlo_process;

    unsigned int host_nparts_;

private:
};

#endif
