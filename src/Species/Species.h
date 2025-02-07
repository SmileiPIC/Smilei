#ifndef SPECIES_H
#define SPECIES_H

#include <vector>
#include <string>
// #include "PyTools.h"

#include "Particles.h"
#ifdef SMILEI_ACCELERATOR_GPU_OACC
#include "nvidiaParticles.h"
#endif
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
#include "Merging.h"
#include "PartCompTime.h"
#include "BirthRecords.h"

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
class Merging;
class PartCompTime;


//! class Species
class Species
{
public:
    // -----------------------------------------------------------------------------
    //  1. Constructor/destructor

    //! Species constructor
    Species( Params &, Patch * );

    //! Species destructor
    virtual ~Species();

    // -----------------------------------------------------------------------------
    //  2. Species parameters

    //! number of the species
    unsigned int species_number_;

    //! kind/name of species
    std::string name_;

    //! position initialization type, possible values: "regular" or "random"
    std::string position_initialization_;

    //! momentum initialization type, possible values: "cold" or "maxwell-juettner"
    std::string momentum_initialization_;

    //! coefficient on the maximum number of particles for the species
    double c_part_max_;

    //! mass [electron mass]
    double mass_;

    //! atomic number
    unsigned int atomic_number_;

    //! maximum charge state
    unsigned int maximum_charge_state_;

    //! alpha parameter in the Tong-Lin ionization model
    double ionization_tl_parameter_;

    //! user defined ionization rate profile
    PyObject *ionization_rate_;

    //! thermalizing temperature for thermalizing BCs [\f$m_e c^2\f$]
    std::vector<double> thermal_boundary_temperature_;
    //! mean velocity used when thermalizing BCs are used [\f$c\f$]
    std::vector<double> thermal_boundary_velocity_;

    //! thermal velocity [\f$c\f$]
    std::vector<double> thermal_velocity_;
    //! thermal momentum [\f$m_e c\f$]
    std::vector<double> thermal_momentum_;

    //! pusher name
    std::string pusher_name_;

    //! radiation model
    std::string radiation_model_;

    //! Time for which the species is frozen
    double time_frozen_;

    //! logical true if particles radiate
    bool radiating_;

    //! logical true if particles are relativistic and require proper electromagnetic field initialization
    bool relativistic_field_initialization_;

    //! Iteration for which the species field is initialized in case of relativistic initialization
    int iter_relativistic_initialization_;

    //! Boundary conditions for particules
    std::vector<std::vector<std::string> > boundary_conditions_;

    //! Ionization model per Specie (tunnel)
    std::string ionization_model_;

    //! Type of density profile ("nb" or "charge")
    std::string density_profile_type_;

    //! charge profile
    Profile *charge_profile_;

    //! density profile
    Profile *density_profile_;

    //! vector of velocity profiles (vx, vy, vz)
    std::vector<Profile *> velocity_profile_;
    
    //! True if velocity profile is radial instead of cartesian
    bool radial_velocity_profile_;

    //! vector of temperature profiles (Tx, Ty, Tz)
    std::vector<Profile *> temperature_profile_;

    //! number-of-particle-per-cell profile
    Profile *particles_per_cell_profile_;

    // -----------------------------------------------------------------------------
    //  3. Variables for species processing

    SpeciesMPIbuffers MPI_buffer_;

    //! Maximum charge at initialization
    double max_charge_;

    //! Vector containing all Particles of the considered Species
    Particles *particles;
    Particles particles_sorted[2];
    //std::vector<int> index_of_particles_to_exchange;

    //! If initialization from file, this contains the number of particles. Otherwise 0
    unsigned int file_position_npart_;
    unsigned int file_momentum_npart_;

    //! Pointer toward position array
    void *position_initialization_array_;
    //! Pointer toward regular number of particles array
    std::vector<int> regular_number_array_;
    //! Number of particles in the init array
    void *momentum_initialization_array_;
    //! Number of particles in the init array
    unsigned int n_numpy_particles_;
    //! Boolean to know if we initialize particles one specie on another species
    bool position_initialization_on_species_;
    //! Index of the species where position initialization is made
    int position_initialization_on_species_index_;
    //! Pointer to the species where field-ionized electrons go
    Species *electron_species;
    //! Index of the species where field-ionized electrons go
    size_t electron_species_index;
    //! Name of the species where field-ionized electrons go
    std::string ionization_electrons;

    //! Pointer to the species where radiated photon go
    Species * photon_species_;
    //! Index of the species where radiated photons go
    int photon_species_index;
    //! radiation photon species for the Monte-Carlo model.
    //! Name of the species where radiated photons go
    std::string radiation_photon_species;
    //! Number of photons emitted per particle and per event
    int radiation_photon_sampling_;
    //! Maximum number of photon emissions process during a timestep per particle
    int radiation_max_emissions_;
    //! Threshold on the photon Lorentz factor under which the macro-photon
    //! is not generated but directly added to the energy scalar diags
    //! This enable to limit emission of useless low-energy photons
    double radiation_photon_gamma_threshold_;
    //! Particles object to store emitted photons by radiation at each time step
    Particles * radiated_photons_ = nullptr;

    //! Pointer to the species where electron-positron pairs
    //! from the multiphoton Breit-Wheeler go
    Species * mBW_pair_species_[2] = {nullptr, nullptr};
    //! Index of the species where electron-positron pairs
    //! from the multiphoton Breit-Wheeler go
    int mBW_pair_species_index_[2];
    //! Number of created pairs per event and per photons
    int mBW_pair_creation_sampling_[2];
    //! electron and positron Species for the multiphoton Breit-Wheeler
    std::vector<std::string> mBW_pair_species_names_;
    // Particles object to store created electron-positron pairs
    Particles * mBW_pair_particles_[2] = {nullptr , nullptr};

    //! Cluster width in number of cells
    unsigned int cluster_width_; //Should divide the number of cells in X of a single MPI domain.
    //! Array counting the occurence of each cell key
    std::vector<int> count;
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

    //! Inverse of the number of spatial dimension for the fields
    double inv_nDim_field;

    //! sub primal dimensions of fields
    unsigned int f_dim0, f_dim1, f_dim2;

    //! sub dual dimensions of fields
    unsigned int f_dim0_d, f_dim1_d, f_dim2_d;

    //! Accumulate energy lost with bc
    double nrj_bc_lost;
    //! Accumulate energy lost with moving window
    double nrj_mw_out;
    //! Accumulate energy gained with moving window
    double nrj_mw_inj;
    //! Accumulate energy added with new particles
    double nrj_new_part_;
    //! Accumulate energy lost by the particle with the radiation
    double nrj_radiated_;

    //! whether to choose vectorized operators with respective sorting methods
    int vectorized_operators;

    // Merging parameters :
    //! Merging method
    std::string merging_method_;

    //! Boolean to test if the species has the merging ready
    bool has_merging_;

    //! Time selection for the particle merging
    TimeSelection *merging_time_selection_;

    //! Minimum number of particles in a packet that can be merged
    //! at the same time
    unsigned int merge_min_packet_size_;

    //! Maximum number of particles in a packet that can be merged
    //! at the same time
    unsigned int merge_max_packet_size_;

    //! Minimum number of particles per cell to be able to merge
    unsigned int merge_min_particles_per_cell_;

    //! Flag to activate the correction against the accumulation effect
    bool merge_accumulation_correction_;

    //! Minimum momentum cell length for the merging
    std::vector<double>  merge_min_momentum_cell_length_;

    //! Discreatization in the momentum space (number of momentum cells in each direction)
    std::vector<unsigned int> merge_momentum_cell_size_;

    //! Discretization scale
    bool merge_log_scale_;

    //! Minimum momentum value in log scale
    double merge_min_momentum_log_scale_;

    //! Local minimum of MPI domain
    double min_loc;

    //! Inverse of the number of spatial dimension for the particles
    double inv_nDim_particles;

    // -----------------------------------------------------------------------------
    //  4. Operators

    //! Ionization method
    Ionization *Ionize;

    //! Radiation method (Continuous or Monte-Carlo)
    Radiation *Radiate;

    //! Multiphoton Breit-wheeler
    MultiphotonBreitWheeler *Multiphoton_Breit_Wheeler_process;

    //! Boundary condition for the Particles of the considered Species
    PartBoundCond *partBoundCond;

    //! Particles pusher (change momentum & change position, only momentum in case envelope model is used)
    Pusher *Push;

    //! Particles position pusher (change change position)
    Pusher *Push_ponderomotive_position = NULL;

    //! Interpolator (used to push particles and for probes)
    Interpolator *Interp;

    //! Projector
    Projector *Proj;

    //! Merging
    Merging *Merge;

    //! Particle Computation time evaluation
    PartCompTime *part_comp_time_ = NULL;
    
    //! Birth records
    BirthRecords *birth_records_ = NULL;

    // -----------------------------------------------------------------------------
    //  5. Methods

    virtual void initCluster( Params &, Patch * );

    virtual void resizeCluster( Params & );

    //! Initialize particles
    void initParticles( Params &params, Patch *patch, bool with_particles = true, Particles * like_particles = NULL );

    //! Initialize operators (must be separate from parameters init, because of cloning)
    void initOperators( Params &, Patch * );

    //! Method returning the Particle list for the considered Species
    inline Particles getParticlesList() const
    {
        return *particles;
    }

    //! Method returning the Particle list pointer for the considered Species
    inline Particles &getParticlesList()
    {
        return *particles;
    }

    //! Method returning the effective number of Particles for the considered Species
    inline unsigned int getNbrOfParticles() const
    {
        return particles->numberOfParticles();
    }

    //! Method returning the size of Particles
    inline unsigned int getParticlesSize() const
    {
        return particles->size();
    }

    //! Return the capacity of Particles
    inline unsigned int getParticlesCapacity() const
    {
        return particles->capacity();
    }

#if defined( SMILEI_ACCELERATOR_GPU )

    void allocateParticlesOnDevice();

    //! Copy particles from host to device
    void
    copyParticlesFromHostToDevice();

    //! Prepare the species Current and Rho grids on Device
    void
    prepareSpeciesCurrentAndChargeOnDevice( 
        unsigned int ispec,
        ElectroMagn * EMfields
    );

    //! Deallocate species Current (J) and Charge (Rho) arrays on Device
    void
    deleteSpeciesCurrentAndChargeOnDevice(
        unsigned int ispec,
        ElectroMagn * EMfields
    );

#endif

    //! Method calculating the Particle dynamics (interpolation, pusher, projection and more)
    //! For all particles of the species
    //!   - interpolate the fields at the particle position
    //!   - perform ionization
    //!   - perform the radiation reaction
    //!   - calculate the new velocity
    //!   - calculate the new position
    //!   - apply the boundary conditions
    //!   - increment the currents (projection)
    virtual void dynamics( double time, unsigned int ispec,
                           ElectroMagn *EMfields,
                           Params &params, bool diag_flag,
                           PartWalls *partWalls, Patch *patch, SmileiMPI *smpi,
                           RadiationTables &RadiationTables,
                           MultiphotonBreitWheelerTables &MultiphotonBreitWheelerTables );

    //! Method projecting susceptibility and calculating the particles updated momentum (interpolation, momentum pusher), only particles interacting with envelope
    virtual void ponderomotiveUpdateSusceptibilityAndMomentum( double time_dual,
            ElectroMagn *EMfields,
            Params &params,
            Patch *patch, SmileiMPI *smpi );
    //! Method projecting susceptibility, only particles interacting with envelope
    virtual void ponderomotiveProjectSusceptibility( double time_dual,
            ElectroMagn *EMfields,
            Params &params,
            Patch *patch, SmileiMPI *smpi );

    //! Method calculating the Particle updated position (interpolation, position pusher, only particles interacting with envelope)
    // and projecting charge density and thus current density (through Esirkepov method) for Maxwell's Equations
    virtual void ponderomotiveUpdatePositionAndCurrents( double time_dual, unsigned int ispec,
            ElectroMagn *EMfields,
            Params &params, bool diag_flag, PartWalls *partWalls,
            Patch *patch, SmileiMPI *smpi );

    //! Method calculating the Particle dynamics with scalar operators (interpolation, pusher, projection)
    virtual void scalarDynamics( double time, unsigned int ispec,
                                  ElectroMagn *EMfields,
                                  Params &params, bool diag_flag,
                                  PartWalls *partWalls, Patch *patch, SmileiMPI *smpi,
                                  RadiationTables &RadiationTables,
                                  MultiphotonBreitWheelerTables &MultiphotonBreitWheelerTables );

    virtual void scalarPonderomotiveUpdateSusceptibilityAndMomentum( double, 
            ElectroMagn *,
            Params &, 
            Patch *, SmileiMPI * ) {};

    virtual void scalarPonderomotiveUpdatePositionAndCurrents( double, unsigned int,
            ElectroMagn *,
            Params &, bool, PartWalls *,
            Patch *, SmileiMPI * ) {};

    //! Projection method used specifically for the diagnotics
    virtual void projectionForDiags(  unsigned int ispec,
                                       ElectroMagn *EMfields,
                                       Params &params, bool diag_flag,
                                       Patch *patch );

    //! Method performing the importation of new particles
    virtual void dynamicsImportParticles( double time, 
                                            Params &params,
                                            Patch *patch,
                                            std::vector<Diagnostic *> &localDiags );

    //! Method performing the merging of particles
    virtual void mergeParticles( double time_dual );


    //! Method calculating the Particle charge on the grid (projection)
    virtual void computeCharge( ElectroMagn *EMfields, bool old=false );

    //! Method used to inject and sort particles
    virtual void sortParticles( Params &param );

    virtual void computeParticleCellKeys(   Params    &,
                                            Particles *,
                                            int       * __restrict__,
                                            int       * __restrict__,
                                            unsigned int,
                                            unsigned int ) {};

    virtual void computeParticleCellKeys( Params & ) {};

    //! This function configures the type of species according to the default mode
    //! regardless the number of particles per cell
    virtual void defaultConfigure( Params &params, Patch *patch ) ;

    //! This function configures the species according to the vectorization mode
    virtual void configuration( Params &params, Patch *patch ) ;

    //! This function reconfigures the species operators after evaluating
    //! the best mode from the particle distribution
    virtual void reconfiguration( Params &param, Patch   *patch );

    //! Counting sort method for particles
    void countSortParticles( Params &param );

    //!
    virtual void addSpaceForOneParticle()
    {
        particles->last_index[particles->last_index.size()-1]++;
    }

    //! Method to know if we have to project this species or not.
    bool  isProj( double time_dual, SimWindow *simWindow );

    inline double computeEnergy()
    {
        double nrj( 0. );
        if( mass_ > 0 ) {
            for( unsigned int iPart=0 ; iPart<getNbrOfParticles() ; iPart++ ) {
                nrj += particles->weight( iPart )*( particles->LorentzFactor( iPart )-1.0 );
            }
            nrj *= mass_;
        } else if( mass_ == 0 ) {
            for( unsigned int iPart=0 ; iPart<getNbrOfParticles() ; iPart++ ) {
                nrj += particles->weight( iPart )*( particles->momentumNorm( iPart ) );
            }
        }
        return nrj;
    }

    inline std::size_t getMemFootPrint()
    {
        /*int speciesSize  = ( 2*nDim_particle + 3 + 1 )*sizeof(double) + sizeof(short);
        if ( particles->is_test )
            speciesSize += sizeof ( unsigned int );*/
        //speciesSize *= getNbrOfParticles();
        std::size_t speciesSize = 0;
        speciesSize += particles->double_prop_.size()*sizeof( double );
        speciesSize += particles->short_prop_.size()*sizeof( short );
        speciesSize += particles->uint64_prop_.size()*sizeof( uint64_t );
        speciesSize *= getParticlesCapacity();
        return speciesSize;
    }

    //! Method to import particles in this species while conserving the sorting among bins
    virtual void importParticles( Params &, Patch *, Particles &, std::vector<Diagnostic *> &, double time_dual, Ionization *I = nullptr );

    //! This method eliminates the space gap between the bins
    //! (presence of empty particles between the bins)
    void compress(SmileiMPI *smpi,
        int ithread,
        bool compute_cell_keys = false);

    //! This method removes particles with a negative weight
    //! without changing the bin first index
    //! Bins are therefore potentially seperated by empty particle slots
    void removeTaggedParticlesPerBin(
        SmileiMPI *smpi,
        int ithread,
        bool compute_cell_keys = false);

    //! This method removes particles with a negative weight
    //! when a single bin is used
#ifdef SMILEI_ACCELERATOR_GPU_OACC
    void removeTaggedParticles(
        SmileiMPI *smpi,
        int *const first_index,
        int *const last_index,
        int ithread,
        bool compute_cell_keys = false);
#endif

    //! Moving window boundary conditions managment
    void disableXmax();
    //! Moving window boundary conditions managment
    void setXminBoundaryCondition();

    //! Check function that enables to control the results of some operators
    void check( Patch *patch, std::string title );

    //! Perform the sum of all Lorentz factor
    double sumGamma()
    {
        double s_gamma( 0. );
        for( unsigned int ipart = 0 ; ipart < getNbrOfParticles() ; ipart++ ) {
            s_gamma += sqrt( 1. + particles->momentum( 0, ipart ) * particles->momentum( 0, ipart )
                             + particles->momentum( 1, ipart ) * particles->momentum( 1, ipart )
                             + particles->momentum( 2, ipart ) * particles->momentum( 2, ipart ) );
        }

        return s_gamma;
    }

    typedef double (Species::*fptr)(Particles*, int, int);
    //Array of pointers to functions measuring distance along each dimension.
    fptr distance[3];

    double cartesian_distance(Particles *part, int idim, int ipart){
        return part->position(idim, ipart) - min_loc_vec[idim];
    }

    double radial_distance(Particles *part, int idim, int ipart){
        return sqrt(  part->position(idim  , ipart) * part->position(idim  , ipart)
                    + part->position(idim+1, ipart) * part->position(idim+1, ipart))
               - min_loc_vec[idim];
    }

    //! Erase all particles with zero weight
    void eraseWeightlessParticles();

    // ---- Variables for tasks -> does not look like it's only used for tasks now ...

    // Number of bins for the use of tasks
    unsigned int Nbins;

    // Number of cells used fot tasks + vectorization
    int Ncells;

    // buffers for bin projection when tasks are used
    std::vector<double *> b_Jx;
    std::vector<double *> b_Jy;
    std::vector<double *> b_Jz;
    std::vector<double *> b_rho;
    std::vector<double *> b_Chi;

    // Tags for the task dependencies of the particle operations
    int *bin_has_interpolated;
    int *bin_has_pushed;
    int *bin_has_done_particles_BC;
    int *bin_has_projected_chi;
    int *bin_has_projected;

    // buffers for bin projection when tasks are used
    std::vector< std::complex<double> *> b_Jl;
    std::vector< std::complex<double> *> b_Jr;
    std::vector< std::complex<double> *> b_Jt;
    std::vector< std::complex<double> *> b_rhoAM;
    std::vector<double *> b_ChiAM;

    //! Size of the projection buffer
    unsigned int size_proj_buffer_Jx,size_proj_buffer_Jy,size_proj_buffer_Jz,size_proj_buffer_rho;
    unsigned int size_proj_buffer_Jl,size_proj_buffer_Jr,size_proj_buffer_Jt,size_proj_buffer_rhoAM;
    std::string geometry;
    double *nrj_lost_per_bin;
    double *radiated_energy_per_bin;

protected:

    //! Patch length
    unsigned int length_[3];

private:

};

#endif
