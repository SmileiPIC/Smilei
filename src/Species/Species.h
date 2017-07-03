#ifndef SPECIES_H
#define SPECIES_H

#include <vector>
#include <string>

#include "Particles.h"
#include "Params.h"
//#include "PartBoundCond.h"

#include "Pusher.h"
#include "Ionization.h"
#include "ElectroMagn.h"
#include "Profile.h"
#include "AsyncMPIbuffers.h"

class ElectroMagn;
class Pusher;
class Interpolator;
class Projector;
class PartBoundCond;
class PartWalls;
class Field3D;
class Patch;
class SimWindow;


//! class Species
class Species
{
public:
    SpeciesMPIbuffers MPIbuff;
    
    //! Species creator
    Species(Params&, Patch*);
    
    void initCluster(Params&);
    
    //! Initialize operators (must be separate from parameters init, because of cloning)
    void initOperators(Params&, Patch*);
    
    //! Species destructor
    virtual ~Species();
    
    //! number of the species
    unsigned int speciesNumber;
    
    //! kind/name of species
    std::string name;
    
    //! position initialization type, possible values: "regular" or "random"
    std::string position_initialization;
    
    //! momentum initialization type, possible values: "cold" or "maxwell-juettner"
    std::string momentum_initialization;
    
    //! coefficient on the maximum number of particles for the species
    double c_part_max;
    
    //! mass [electron mass]
    double mass;
    
    //! atomic number
    unsigned int atomic_number;
    
    //! thermalizing temperature for thermalizing BCs [\f$m_e c^2\f$]
    std::vector<double> thermal_boundary_temperature;
    //! mean velocity used when thermalizing BCs are used [\f$c\f$]
    std::vector<double> thermal_boundary_velocity;
    
    //! thermal velocity [\f$c\f$]
    std::vector<double> thermalVelocity;
    //! thermal momentum [\f$m_e c\f$]
    std::vector<double> thermalMomentum;
    
    //! pusher name
    std::string pusher;
    
    //! Time for which the species is frozen
    double time_frozen;
    
    //! logical true if particles radiate
    bool radiating;
    
    //! Boundary conditions for particules
    std::vector<std::vector<std::string> > boundary_conditions;
    
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
        return *particles;
    }
    inline Particles& getParticlesList() {
        return *particles;
    }
    
    //! Method returning the effective number of Particles for the considered Species
    inline unsigned int getNbrOfParticles() const {
        return (*particles).size();
    }
    // capacity() = vect ever oversize
    //! \todo define particles.capacity = min.capacity
    inline unsigned int getParticlesCapacity() const {
        return (*particles).capacity();
    }
    
    //! Method calculating the Particle dynamics (interpolation, pusher, projection)
    virtual void dynamics(double time, unsigned int ispec, ElectroMagn* EMfields, Interpolator* interp,
                          Projector* proj, Params &params, bool diag_flag,
                          PartWalls* partWalls, Patch* patch, SmileiMPI* smpi, std::vector<Diagnostic*>& localDiags);
    
    //! Method calculating the Particle charge on the grid (projection)
    virtual void computeCharge(unsigned int ispec, ElectroMagn* EMfields, Projector* Proj);
    
    //! Method used to initialize the Particle position in a given cell
    void initPosition(unsigned int, unsigned int, double *);
    
    //! Method used to initialize the Particle 3d momentum in a given cell
    void initMomentum(unsigned int, unsigned int, double *, double *);
    
    //! Method used to initialize the Particle weight (equivalent to a charge density) in a given cell
    void initWeight(unsigned int,  unsigned int, double);
    
    //! Method used to initialize the Particle charge
    void initCharge(unsigned int, unsigned int, double);
    
    //! Maximum charge at initialization
    double max_charge;
    
    //! Method used to sort particles
    void sort_part();
    void count_sort_part(Params& param);
    
    //void updateMvWinLimits(double x_moved);
    
    //! Vector containing all Particles of the considered Species
    Particles *particles;
    Particles particles_sorted[2];
    //std::vector<int> index_of_particles_to_exchange;
    
    //! Ionization method
    Ionization* Ionize;
    
    //! Pointer to the species where field-ionized electrons go
    Species *electron_species;
    //! Index of the species where field-ionized electrons go
    int electron_species_index;
    //! Name of the species where field-ionized electrons go
    std::string ionization_electrons;
    
    //! Cluster width in number of cells
    unsigned int clrw; //Should divide the number of cells in X of a single MPI domain. 
    //! first and last index of each particle bin
    std::vector<int> bmin, bmax;
    //! sub dimensions of buffers for dim > 1
    std::vector<unsigned int> b_dim;
    
    //! Oversize (copy from Params)
    std::vector<unsigned int> oversize;
    
    //! MPI structure to exchange particles
    std::vector<MPI_Datatype> typePartSend ;
    std::vector<MPI_Datatype> typePartRecv ;
    MPI_Datatype exchangePatch;

    //! Cell_length (copy from Params)
    std::vector<double> cell_length;
    //! min_loc_vec (copy from picparams)
    std::vector<double> min_loc_vec;
    
    //inline void clearExchList(int tid) {
    //        indexes_of_particles_to_exchange_per_thd[tid].clear();
    //}
    inline void clearExchList() {
        indexes_of_particles_to_exchange.clear();
    }
    //inline void addPartInExchList(int tid, int iPart) {
    //    indexes_of_particles_to_exchange_per_thd[tid].push_back(iPart);
    //}
    inline void addPartInExchList( int iPart) {
        indexes_of_particles_to_exchange.push_back(iPart);
    }
    //std::vector< std::vector<int> > indexes_of_particles_to_exchange_per_thd;
    std::vector<int>                indexes_of_particles_to_exchange;
    //std::vector<int>                new_indexes_of_particles_to_exchange;
    
    //! Method to know if we have to project this species or not.
    bool  isProj(double time_dual, SimWindow* simWindow);
    
    double getLostNrjBC() const {return mass*nrj_bc_lost;}
    double getLostNrjMW() const {return mass*nrj_mw_lost;}
    
    double getNewParticlesNRJ() const {return mass*nrj_new_particles;}
    void reinitDiags() { 
        //nrj_bc_lost = 0;
        nrj_mw_lost = 0;
        nrj_new_particles = 0;
    }
    inline void storeNRJlost( double nrj ) { nrj_mw_lost += nrj; };
    
    inline double computeNRJ() {
        double nrj(0.);
        for ( unsigned int iPart=0 ; iPart<getNbrOfParticles() ; iPart++ )
            nrj += (*particles).weight(iPart)*((*particles).lor_fac(iPart)-1.0);
        return nrj;
    }
    
    inline int getMemFootPrint() {
        /*int speciesSize  = ( 2*nDim_particle + 3 + 1 )*sizeof(double) + sizeof(short);
        if ( particles->is_test )
            speciesSize += sizeof ( unsigned int );*/
        //speciesSize *= getNbrOfParticles();
        int speciesSize(0);
        speciesSize += particles->double_prop.size()*sizeof(double);
        speciesSize += particles->short_prop.size()*sizeof(short);
        speciesSize += particles->uint64_prop.size()*sizeof(uint64_t);
        speciesSize *= getParticlesCapacity();
        return speciesSize;
    }
    
    //! Method to create new particles.
    int  createParticles(std::vector<unsigned int> n_space_to_create, Params& params, Patch * patch, int new_bin_idx);
    
    //! Method to import particles in this species while conserving the sorting among bins
    void importParticles( Params&, Patch*, Particles&, std::vector<Diagnostic*>& );
    
    //! 2 times pi
    double PI2;
    double PI_ov_2;
    double dx_inv_, dy_inv_, dz_inv_;
    
    //! Boundary condition for the Particles of the considered Species
    PartBoundCond* partBoundCond;
    
    //! Particles pusher (change momentum & change position)
    Pusher* Push;
    
    //! Moving window boundary conditions managment
    void disableXmax();
    //! Moving window boundary conditions managment
    void setXminBoundaryCondition();
    
    //! Number of the associated tracking diagnostic
    unsigned int tracking_diagnostic;
    
private:
    //! Number of steps for Maxwell-Juettner cumulative function integration
    //! \todo{Put in a code constant class}
//    unsigned int nE;
    
    //! Parameter used when defining the Maxwell-Juettner function (corresponds to a Maximum energy)
//    double muEmax;
    
    //! Parameter used when defining the Maxwell-Juettner function (corresponds to a discretization step in energy)
//    double dE;
    
    //! Number of spatial dimension for the particles
    unsigned int nDim_particle;
    
    //! Number of spatial dimension for the fields
    unsigned int nDim_field;
    
    //! Inverse of the number of spatial dimension for the fields
    double inv_nDim_field;
    
    //! Local minimum of MPI domain
    double min_loc;
    
    //! sub primal dimensions of fields
    unsigned int f_dim0, f_dim1, f_dim2;
    
    //! Accumulate nrj lost with bc
    double nrj_bc_lost;
    //! Accumulate nrj lost with moving window
    double nrj_mw_lost;
    //! Accumulate nrj added with new particles
    double nrj_new_particles;
    
    //! Samples npoints values of energies in a Maxwell-Juttner distribution
    std::vector<double> maxwellJuttner(unsigned int npoints, double temperature);
    //! Array used in the Maxwell-Juttner sampling (see doc)
    static const double lnInvF[1000];
    //! Array used in the Maxwell-Juttner sampling (see doc)
    static const double lnInvH[1000];

};

#endif
