#ifndef SPECIES_H
#define SPECIES_H

#include <vector>
#include <string>

#include "Particles.h"
#include "Params.h"
//#include "Pusher.h"
//#include "PartBoundCond.h"
#include "Params.h"
#include "SmileiMPI.h"
#include "Pusher.h"
#include "Ionization.h"
#include "ElectroMagn.h"
#include "Profile.h"

class ElectroMagn;
class Pusher;
class Interpolator;
class Projector;
class PartBoundCond;
class PartWall;
class Field3D;
class SimWindow;


//! class Species
class Species
{
public:
    //! Species creator
    Species(Params&, SmileiMPI*);
    //! Species destructor
    virtual ~Species();
    
    //! number of the species
    unsigned int speciesNumber;
    
    //! kind/name of species
    std::string species_type;
    
    //! position initialization type, possible values: "regular" or "random"
    std::string initPosition_type;
    
    //! momentum initialization type, possible values: "cold" or "maxwell-juettner"
    std::string initMomentum_type;
    
    //! coefficient on the maximum number of particles for the species
    double c_part_max;
    
    //! mass [electron mass]
    double mass;
    
    //! atomic number
    unsigned int atomic_number;
    
    //! thermalizing temperature for thermalizing BCs [\f$m_e c^2\f$]
    std::vector<double> thermT;
    //! mean velocity used when thermalizing BCs are used [\f$c\f$]
    std::vector<double> thermVelocity;
    
    //! thermal velocity [\f$c\f$]
    std::vector<double> thermalVelocity;
    //! thermal momentum [\f$m_e c\f$]
    std::vector<double> thermalMomentum;
    
    //! dynamics type. Possible values: "Norm" "Radiation Reaction"
    std::string dynamics_type;
    
    //! Time for which the species is frozen
    double time_frozen;
    
    //! logical true if particles radiate
    bool radiating;
    
    //! logical true if particles radiate
    bool isTest;
    
    //! dump every for test particles
    unsigned int test_dump_every;
    
    //! Boundary conditions for particules
    std::string bc_part_type_west;
    std::string bc_part_type_east;
    std::string bc_part_type_south;
    std::string bc_part_type_north;
    std::string bc_part_type_bottom;
    std::string bc_part_type_up;
    
    //! Ionization model per Specie (tunnel)
    std::string ionization_model;
    
    //! Type of density profile ("nb" or "charge")
    std::string densityProfileType;
    
    //! charge profile
    Profile *chargeProfile;
    
    //! density profile
    Profile *densityProfile;
    
    //! vector of velocity profiles (vx, vy, vz)
    std::vector<Profile *> velocityProfile;
    
    //! vector of temperature profiles (Tx, Ty, Tz)
    std::vector<Profile *> temperatureProfile;
    
    //! number-of-particle-per-cell profile
    Profile *ppcProfile;
    
    
    //! Method returning the Particle list for the considered Species
    inline Particles getParticlesList() const {
        return particles;
    }
    inline Particles& getParticlesList() {
        return particles;
    }
    
    //! Method returning the effective number of Particles for the considered Species
    inline unsigned int getNbrOfParticles() const {
        return particles.size();
    }
    // capacity() = vect ever oversize
    //! \todo define particles.capacity = min.capacity
    inline unsigned int getParticlesCapacity() const {
        return particles.capacity();
    }
    
    //! Method calculating the Particle dynamics (interpolation, pusher, projection)
    virtual void dynamics(double time, ElectroMagn* EMfields, Interpolator* interp,
                          Projector* proj, SmileiMPI *smpi, Params &params, SimWindow* simWindow,
                          std::vector<PartWall*> vecPartWall);
    
    //! Method used to initialize the Particle position in a given cell
    void initPosition(unsigned int, unsigned int, double *, unsigned int);
    
    //! Method used to initialize the Particle 3d momentum in a given cell
    void initMomentum(unsigned int, unsigned int, double *, double *, std::vector<double>&);
    
    //! Method used to initialize the Particle weight (equivalent to a charge density) in a given cell
    void initWeight(unsigned int,  unsigned int, double);
    
    //! Method used to initialize the Particle charge
    void initCharge(unsigned int, unsigned int, double);
    
    //! Maximum charge at initialization
    double max_charge;
    
    //! Method used to save all Particles properties for the considered Species
    void dump(std::ofstream&);
    
    //! Method used to sort particles
    void sort_part();
    
    void movingWindow_x(unsigned int shift, SmileiMPI *smpi, Params& param);
    void defineNewCells(unsigned int shift, SmileiMPI *smpi, Params& param);
    void updateMvWinLimits(double x_moved);
    
    //! Object containing all Particles of the considered Species
    Particles particles;
    //std::vector<int> index_of_particles_to_exchange;
    
    //! Ionization method
    Ionization* Ionize;
    
    //! to keep rack of ionized electrons
    Species *electron_species;
    
    //! Cluster width in number of cells
    unsigned int clrw; //Should divide the number of cells in X of a single MPI domain. Should default to 1.
    //! first and last index of each particle bin
    std::vector<int> bmin, bmax;
    
    //! Oversize (copy from Params)
    std::vector<unsigned int> oversize;
    
    //! Cell_length (copy from Params)
    std::vector<double> cell_length;
    
    inline void clearExchList(int tid) {
        indexes_of_particles_to_exchange_per_thd[tid].clear();
    }
    inline void addPartInExchList(int tid, int iPart) {
        indexes_of_particles_to_exchange_per_thd[tid].push_back(iPart);
    }
    std::vector< std::vector<int> > indexes_of_particles_to_exchange_per_thd;
    
    //! Method to know if we have to project this species or not.
    bool  isProj(double time_dual, SimWindow* simWindow);
    
    double getLostNrjBC() const {return mass*nrj_bc_lost;}
    double getLostNrjMW() const {return mass*nrj_mw_lost;}
    
    double getNewParticlesNRJ() const {return mass*nrj_new_particles;}
    void reinitDiags() { 
    nrj_bc_lost = 0;
    nrj_mw_lost = 0;
    nrj_new_particles = 0;
    }
    
    inline int getMemFootPrint() {
        int speciesSize  = ( 2*ndim + 3 + 1 )*sizeof(double) + sizeof(short);
        if ( particles.isTestParticles )
            speciesSize += sizeof ( unsigned int );
        //speciesSize *= getNbrOfParticles();
        speciesSize *= getParticlesCapacity();
        return speciesSize;
    }
    
    //! Method to create new particles.
    int  createParticles(std::vector<unsigned int> n_space_to_create, std::vector<double> cell_index, int new_bin_idx,  Params& param);
    
    //! Boundary condition for the Particles of the considered Species
    PartBoundCond* partBoundCond;
    
    //! Particles pusher (change momentum & change position)
    Pusher* Push;

private:    
    //! Number of steps for Maxwell-Juettner cumulative function integration
    //! \todo{Put in a code constant class}
    unsigned int nE;
    
    //! Parameter used when defining the Maxwell-Juettner function (corresponds to a Maximum energy)
    double muEmax;
    
    //! Parameter used when defining the Maxwell-Juettner function (corresponds to a discretization step in energy)
    double dE;
    
    //! Number of spatial dimension for the particles
    unsigned int ndim;
    
    //! Local minimum of MPI domain
    double min_loc;
    
    //! Size of the projection buffer
    unsigned int size_proj_buffer;
    
    //! sub dimensions of buffers for dim > 1
    unsigned int b_dim0, b_dim1, b_dim2, b_lastdim;
    //! sub primal dimensions of fields
    unsigned int f_dim0, f_dim1, f_dim2;
    
    //! Accumulate nrj lost with bc
    double nrj_bc_lost;
    //! Accumulate nrj lost with moving window
    double nrj_mw_lost;
    //! Accumulate nrj added with new particles
    double nrj_new_particles;
};

#endif
