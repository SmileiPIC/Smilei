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
#include "Radiation.h"
#include "RadiationTables.h"
#include "MultiphotonBreitWheeler.h"
#include "MultiphotonBreitWheelerTables.h"

class ElectroMagn;
class Pusher;
class Interpolator;
class Projector;
class PartBoundCond;
class PartWalls;
class Field3D;
class Patch;
class SimWindow;
class Radiation;


//! class Species
class Species
{
public:
    // -----------------------------------------------------------------------------
    //  1. Constructor/destructor

    //! Species constructor
    Species(Params&, Patch*);

    //! Species destructor
    virtual ~Species();

    // -----------------------------------------------------------------------------
    //  2. Species parameters

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

    //! radiation model
    std::string radiation_model;

    //! Time for which the species is frozen
    double time_frozen;

    //! logical true if particles radiate
    bool radiating;

    //! logical true if particles are relativistic and require proper electromagnetic field initialization
    bool relativistic_field_initialization;

    //! Time for which the species field is initialized in case of relativistic initialization
    double time_relativistic_initialization;

    //! electron and positron Species for the multiphoton Breit-Wheeler
    std::vector<std::string> multiphoton_Breit_Wheeler;

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

    // -----------------------------------------------------------------------------
    //  3. Variables for species processing

    SpeciesMPIbuffers MPIbuff;

    //! Maximum charge at initialization
    double max_charge;

    //! Vector containing all Particles of the considered Species
    Particles *particles;
    Particles particles_sorted[2];
    //std::vector<int> index_of_particles_to_exchange;

    //! Pointer toward position array
    double *position_initialization_array;
    //! Number of particles in the init array
    double *momentum_initialization_array;
    //! Number of particles in the init array
    unsigned int n_numpy_particles;
    //! Boolean to know if we initialize particles one specie on another species
    bool position_initialization_on_species;
    //! Index of the species where position initialization is made
    int position_initialization_on_species_index;
    //! Boolean to know if species follows ponderomotive loop (laser modeled with envelope)
    bool ponderomotive_dynamics;
    //! Pointer to the species where field-ionized electrons go
    Species *electron_species;
    //! Index of the species where field-ionized electrons go
    int electron_species_index;
    //! Name of the species where field-ionized electrons go
    std::string ionization_electrons;

    //! Pointer to the species where radiated photon go
    Species *photon_species;
    //! Index of the species where radiated photons go
    int photon_species_index;
    //! radiation photon species for the Monte-Carlo model.
    //! Name of the species where radiated photons go
    std::string radiation_photon_species;
    //! Number of photons emitted per particle and per event
    int radiation_photon_sampling;
    //! Threshold on the photon Lorentz factor under which the macro-photon
    //! is not generated but directly added to the energy scalar diags
    //! This enable to limit emission of useless low-energy photons
    double radiation_photon_gamma_threshold;

    //! Pointer to the species where electron-positron pairs
    //! from the multiphoton Breit-Wheeler go
    Species * mBW_pair_species[2];
    //! Index of the species where electron-positron pairs
    //! from the multiphoton Breit-Wheeler go
    int mBW_pair_species_index[2];
    //! Number of created pairs per event and per photons
    std::vector<int> mBW_pair_creation_sampling;

    //! Cluster width in number of cells
    unsigned int clrw; //Should divide the number of cells in X of a single MPI domain.
    //! first and last index of each particle bin
    std::vector<int> bmin, bmax;
    //!
    std::vector<int> species_loc_bmax;
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

    //! 2 times pi
    double PI2;
    double PI_ov_2;
    double dx_inv_[3];

    //! Number of the associated tracking diagnostic
    unsigned int tracking_diagnostic;

    //! Number of spatial dimension for the particles
    unsigned int nDim_particle;

    //! Number of spatial dimension for the fields
    unsigned int nDim_field;

    //! sub primal dimensions of fields
    unsigned int f_dim0, f_dim1, f_dim2;

    //! Accumulate nrj lost with bc
    double nrj_bc_lost;
    //! Accumulate nrj lost with moving window
    double nrj_mw_lost;
    //! Accumulate nrj added with new particles
    double nrj_new_particles;

    //! whether to choose vectorized operators with respective sorting methods
    bool vectorized_operators;

    // -----------------------------------------------------------------------------
    //  4. Operators

    //! Ionization method
    Ionization* Ionize;

    //! Radiation method (Continuous or Monte-Carlo)
    Radiation * Radiate;

    //! Multiphoton Breit-wheeler
    MultiphotonBreitWheeler * Multiphoton_Breit_Wheeler_process;

    //! Boundary condition for the Particles of the considered Species
    PartBoundCond* partBoundCond;

    //! Particles pusher (change momentum & change position, only momentum in case envelope model is used)
    Pusher* Push;

    //! Particles position pusher (change change position)
    Pusher* Push_ponderomotive_position = NULL;


    //! Interpolator (used to push particles and for probes)
    Interpolator* Interp;

    //! Projector
    Projector* Proj;

    // -----------------------------------------------------------------------------
    //  5. Methods

    virtual void initCluster(Params&);

    virtual void resizeCluster(Params&);

    //! Initialize operators (must be separate from parameters init, because of cloning)
    void initOperators(Params&, Patch*);

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
    virtual void dynamics(double time, unsigned int ispec,
                          ElectroMagn* EMfields,
                          Params &params, bool diag_flag,
                          PartWalls* partWalls, Patch* patch, SmileiMPI* smpi,
                          RadiationTables &RadiationTables,
                          MultiphotonBreitWheelerTables & MultiphotonBreitWheelerTables,
                          std::vector<Diagnostic*>& localDiags);

    //! Method projecting susceptibility and calculating the particles updated momentum (interpolation, momentum pusher), only particles interacting with envelope
    virtual void ponderomotive_update_susceptibility_and_momentum(double time_dual, unsigned int ispec,
                           ElectroMagn* EMfields, Interpolator* Interp_envelope,
                           Params &params, bool diag_flag,
                           Patch* patch, SmileiMPI* smpi,
                           std::vector<Diagnostic*>& localDiags);
    //! Method projecting susceptibility, only particles interacting with envelope
    virtual void ponderomotive_project_susceptibility(double time_dual, unsigned int ispec,
                           ElectroMagn* EMfields, Interpolator* Interp_envelope,
                           Params &params, bool diag_flag,
                           Patch* patch, SmileiMPI* smpi,
                           std::vector<Diagnostic*>& localDiags);

    //! Method calculating the Particle updated position (interpolation, position pusher, only particles interacting with envelope)
    // and projecting charge density and thus current density (through Esirkepov method) for Maxwell's Equations
    virtual void ponderomotive_update_position_and_currents(double time_dual, unsigned int ispec,
                           ElectroMagn* EMfields, Interpolator* Interp_envelope, Projector* Proj,
                           Params &params, bool diag_flag, PartWalls* partWalls,
                           Patch* patch, SmileiMPI* smpi,
                           std::vector<Diagnostic*>& localDiags);
    //! Method calculating the Particle dynamics (interpolation, pusher, projection)
    virtual void scalar_dynamics(double time, unsigned int ispec,
                        ElectroMagn* EMfields,
                        Params &params, bool diag_flag,
                        PartWalls* partWalls, Patch* patch, SmileiMPI* smpi,
                        RadiationTables &RadiationTables,
                        MultiphotonBreitWheelerTables & MultiphotonBreitWheelerTables,
                        std::vector<Diagnostic*>& localDiags);

    virtual void projection_for_diags(double time, unsigned int ispec,
                          ElectroMagn* EMfields,
                          Params &params, bool diag_flag,
                          Patch* patch, SmileiMPI* smpi);
    
    //! Method performing the importation of new particles
    virtual void dynamics_import_particles(double time, unsigned int ispec,
                        Params &params,
                        Patch* patch, SmileiMPI* smpi,
                        RadiationTables &RadiationTables,
                        MultiphotonBreitWheelerTables & MultiphotonBreitWheelerTables,
                        std::vector<Diagnostic*>& localDiags);

    //! Method calculating the Particle charge on the grid (projection)
    virtual void computeCharge(unsigned int ispec, ElectroMagn* EMfields);

    //! Method used to initialize the Particle position in a given cell
    void initPosition(unsigned int, unsigned int, double *);

    //! Method used to initialize the Particle 3d momentum in a given cell
    void initMomentum(unsigned int, unsigned int, double *, double *);

    //! Method used to initialize the Particle weight (equivalent to a charge density) in a given cell
    void initWeight(unsigned int,  unsigned int, double);

    //! Method used to initialize the Particle charge
    void initCharge(unsigned int, unsigned int, double);

    //! Method used to sort particles
    virtual void sort_part(Params& param);

    //! This function configures the species according to the vectorization mode
    virtual void configuration( Params& params, Patch * patch) ;

    virtual void reconfiguration(Params& param, Patch  * patch);

    void count_sort_part(Params& param);

    //!
    virtual void add_space_for_a_particle() {
        bmax[bmax.size()-1]++;
    }

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

    //! Get the energy lost in the boundary conditions
    double getLostNrjBC() const {return mass*nrj_bc_lost;}

    //! Get energy lost with moving window (fields)
    double getLostNrjMW() const {return mass*nrj_mw_lost;}

    //! Get the energy radiated away by the particles
    double getNrjRadiation() const {return nrj_radiation;}

    //! Set the energy radiated away by the particles
    void setNrjRadiation(double value) {nrj_radiation = value;}

    //! Add the energy radiated away by the particles
    void addNrjRadiation(double value) {nrj_radiation += value;}

    //! Get energy gained via new particles
    double getNewParticlesNRJ() const {return mass*nrj_new_particles;}

    //! Reinitialize the scalar diagnostics buffer
    void reinitDiags() {
        //nrj_bc_lost = 0;
        nrj_mw_lost = 0;
        nrj_new_particles = 0;
        //nrj_radiation = 0;
    }

    inline void storeNRJlost( double nrj ) { nrj_mw_lost += nrj; };

    inline double computeNRJ() {
        double nrj(0.);
        if (this->mass > 0)
        {
            for ( unsigned int iPart=0 ; iPart<getNbrOfParticles() ; iPart++ )
                nrj += (*particles).weight(iPart)*((*particles).lor_fac(iPart)-1.0);
        }
        else if (this->mass == 0)
        {
            for ( unsigned int iPart=0 ; iPart<getNbrOfParticles() ; iPart++ )
                nrj += (*particles).weight(iPart)*((*particles).momentum_norm(iPart));
        }
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
    virtual void importParticles( Params&, Patch*, Particles&, std::vector<Diagnostic*>& );

    //! Moving window boundary conditions managment
    void disableXmax();
    //! Moving window boundary conditions managment
    void setXminBoundaryCondition();

    //! Check function that enables to control the results of some operators
    void check(Patch * patch, std::string title);

    //! Perform the sum of all Lorentz factor
    double sum_gamma () {
        double s_gamma(0.);
        for ( unsigned int ipart = 0 ; ipart < getNbrOfParticles() ; ipart++ ) {
            s_gamma += sqrt( 1. + particles->momentum(0,ipart) * particles->momentum(0,ipart)
                             + particles->momentum(1,ipart) * particles->momentum(1,ipart)
                             + particles->momentum(2,ipart) * particles->momentum(2,ipart) );
        }

        return s_gamma;
    }


protected:

    //! Accumulate nrj lost by the particle with the radiation
    double nrj_radiation;

    //! Patch length
    unsigned int length[3];

private:
    //! Number of steps for Maxwell-Juettner cumulative function integration
    //! \todo{Put in a code constant class}

//    unsigned int nE;

    //! Parameter used when defining the Maxwell-Juettner function (corresponds to a Maximum energy)
//    double muEmax;

    //! Parameter used when defining the Maxwell-Juettner function (corresponds to a discretization step in energy)
//    double dE;

    //! Inverse of the number of spatial dimension for the fields
    double inv_nDim_field;

    //! Local minimum of MPI domain
    double min_loc;

    //! Samples npoints values of energies in a Maxwell-Juttner distribution
    std::vector<double> maxwellJuttner(unsigned int npoints, double temperature);
    //! Array used in the Maxwell-Juttner sampling (see doc)
    static const double lnInvF[1000];
    //! Array used in the Maxwell-Juttner sampling (see doc)
    static const double lnInvH[1000];

};

#endif
