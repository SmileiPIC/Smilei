#ifndef SPECIES_H
#define SPECIES_H

#include <vector>
#include <string>
//#include "PyTools.h"

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
#include "Merging.h"

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

    //! Time for which the species field is initialized in case of relativistic initialization
    double time_relativistic_initialization_;

    //! electron and positron Species for the multiphoton Breit-Wheeler
    std::vector<std::string> multiphoton_Breit_Wheeler_;

    //! Boundary conditions for particules
    std::vector<std::vector<std::string> > boundary_conditions;

    //! Ionization model per Specie (tunnel)
    std::string ionization_model;

    //! Type of density profile ("nb" or "charge")
    std::string density_profile_type_;

    //! charge profile
    Profile *charge_profile_;

    //! density profile
    Profile *density_profile_;

    //! vector of velocity profiles (vx, vy, vz)
    std::vector<Profile *> velocity_profile_;

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

    //! Pointer toward position array
    double *position_initialization_array_;
    //! Pointer toward regular number of particles array
    std::vector<int> regular_number_array_;
    //! Number of particles in the init array
    double *momentum_initialization_array_;
    //! Number of particles in the init array
    unsigned int n_numpy_particles_;
    //! Boolean to know if we initialize particles one specie on another species
    bool position_initialization_on_species_;
    //! Index of the species where position initialization is made
    int position_initialization_on_species_index;
    //! Initialization type of the species where position initialization is made
    std::string position_initialization_on_species_type_;
    //! Boolean to know if species follows ponderomotive loop (laser modeled with envelope)
    bool ponderomotive_dynamics;
    //! Pointer to the species where field-ionized electrons go
    Species *electron_species;
    //! Index of the species where field-ionized electrons go
    int electron_species_index;
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
    //! Threshold on the photon Lorentz factor under which the macro-photon
    //! is not generated but directly added to the energy scalar diags
    //! This enable to limit emission of useless low-energy photons
    double radiation_photon_gamma_threshold_;

    //! Pointer to the species where electron-positron pairs
    //! from the multiphoton Breit-Wheeler go
    Species *mBW_pair_species[2];
    //! Index of the species where electron-positron pairs
    //! from the multiphoton Breit-Wheeler go
    int mBW_pair_species_index[2];
    //! Number of created pairs per event and per photons
    std::vector<int> mBW_pair_creation_sampling_;

    //! Cluster width in number of cells
    unsigned int clrw; //Should divide the number of cells in X of a single MPI domain.
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

    //! Accumulate energy lost with bc
    double nrj_bc_lost;
    //! Accumulate energy lost with moving window
    double nrj_mw_lost;
    //! Accumulate energy added with new particles
    double new_particles_energy_;
    //! Accumulate energy lost by the particle with the radiation
    double nrj_radiation;

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

    // -----------------------------------------------------------------------------
    //  5. Methods

    virtual void initCluster( Params & );

    virtual void resizeCluster( Params & );

    //! Initialize operators (must be separate from parameters init, because of cloning)
    void initOperators( Params &, Patch * );

    //! Method returning the Particle list for the considered Species
    inline Particles getParticlesList() const
    {
        return *particles;
    }
    inline Particles &getParticlesList()
    {
        return *particles;
    }

    //! Method returning the effective number of Particles for the considered Species
    inline unsigned int getNbrOfParticles() const
    {
        return particles->size();
    }
    // capacity() = vect ever oversize
    //! \todo define particles.capacity = min.capacity
    inline unsigned int getParticlesCapacity() const
    {
        return particles->capacity();
    }

    //! Method calculating the Particle dynamics (interpolation, pusher, projection)
    virtual void dynamics( double time, unsigned int ispec,
                           ElectroMagn *EMfields,
                           Params &params, bool diag_flag,
                           PartWalls *partWalls, Patch *patch, SmileiMPI *smpi,
                           RadiationTables &RadiationTables,
                           MultiphotonBreitWheelerTables &MultiphotonBreitWheelerTables,
                           std::vector<Diagnostic *> &localDiags );

    //! Method projecting susceptibility and calculating the particles updated momentum (interpolation, momentum pusher), only particles interacting with envelope
    virtual void ponderomotiveUpdateSusceptibilityAndMomentum( double time_dual, unsigned int ispec,
            ElectroMagn *EMfields,
            Params &params, bool diag_flag,
            Patch *patch, SmileiMPI *smpi,
            std::vector<Diagnostic *> &localDiags );
    //! Method projecting susceptibility, only particles interacting with envelope
    virtual void ponderomotiveProjectSusceptibility( double time_dual, unsigned int ispec,
            ElectroMagn *EMfields,
            Params &params, bool diag_flag,
            Patch *patch, SmileiMPI *smpi,
            std::vector<Diagnostic *> &localDiags );

    //! Method calculating the Particle updated position (interpolation, position pusher, only particles interacting with envelope)
    // and projecting charge density and thus current density (through Esirkepov method) for Maxwell's Equations
    virtual void ponderomotiveUpdatePositionAndCurrents( double time_dual, unsigned int ispec,
            ElectroMagn *EMfields,
            Params &params, bool diag_flag, PartWalls *partWalls,
            Patch *patch, SmileiMPI *smpi,
            std::vector<Diagnostic *> &localDiags );

    //! Method calculating the Particle dynamics with scalar operators (interpolation, pusher, projection)
    virtual void scalarDynamics( double time, unsigned int ispec,
                                  ElectroMagn *EMfields,
                                  Params &params, bool diag_flag,
                                  PartWalls *partWalls, Patch *patch, SmileiMPI *smpi,
                                  RadiationTables &RadiationTables,
                                  MultiphotonBreitWheelerTables &MultiphotonBreitWheelerTables,
                                  std::vector<Diagnostic *> &localDiags );

    virtual void scalarPonderomotiveUpdateSusceptibilityAndMomentum( double time_dual, unsigned int ispec,
            ElectroMagn *EMfields,
            Params &params, bool diag_flag,
            Patch *patch, SmileiMPI *smpi,
            std::vector<Diagnostic *> &localDiags ) {};

    virtual void scalarPonderomotiveUpdatePositionAndCurrents( double time_dual, unsigned int ispec,
            ElectroMagn *EMfields,
            Params &params, bool diag_flag, PartWalls *partWalls,
            Patch *patch, SmileiMPI *smpi,
            std::vector<Diagnostic *> &localDiags ) {};

    //! Projection method used specifically for the diagnotics
    virtual void projectionForDiags( double time, unsigned int ispec,
                                       ElectroMagn *EMfields,
                                       Params &params, bool diag_flag,
                                       Patch *patch, SmileiMPI *smpi );

    //! Method performing the importation of new particles
    virtual void dynamicsImportParticles( double time, unsigned int ispec,
                                            Params &params,
                                            Patch *patch, SmileiMPI *smpi,
                                            std::vector<Diagnostic *> &localDiags );

    //! Method performing the merging of particles
    virtual void mergeParticles( double time_dual, unsigned int ispec,
                                 Params &params,
                                 Patch *patch, SmileiMPI *smpi,
                                 std::vector<Diagnostic *> &localDiags );


    //! Method calculating the Particle charge on the grid (projection)
    virtual void computeCharge( unsigned int ispec, ElectroMagn *EMfields );

    //! Method used to sort particles
    virtual void sortParticles( Params &param, Patch * patch );

    virtual void computeParticleCellKeys( Params &params ) {};

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

    //inline void clearExchList(int tid) {
    //        indexes_of_particles_to_exchange_per_thd[tid].clear();
    //}
    inline void clearExchList()
    {
        indexes_of_particles_to_exchange.clear();
    }
    //inline void addPartInExchList(int tid, int iPart) {
    //    indexes_of_particles_to_exchange_per_thd[tid].push_back(iPart);
    //}
    inline void addPartInExchList( int iPart )
    {
        indexes_of_particles_to_exchange.push_back( iPart );
    }
    //std::vector< std::vector<int> > indexes_of_particles_to_exchange_per_thd;
    std::vector<int>                indexes_of_particles_to_exchange;
    //std::vector<int>                new_indexes_of_particles_to_exchange;

    //! Method to know if we have to project this species or not.
    bool  isProj( double time_dual, SimWindow *simWindow );

    //! Set the energy lost in the boundary conditions
    void setLostNrjBC( double value )
    {
        nrj_bc_lost = value;
    }
    
    //! Get the energy lost in the boundary conditions
    double getLostNrjBC() const
    {
        return nrj_bc_lost;
    }

    //! Get energy lost with moving window (fields)
    double getLostNrjMW() const
    {
        return nrj_mw_lost;
    }

    //! Get the energy radiated away by the particles
    double getNrjRadiation() const
    {
        return nrj_radiation;
    }

    //! Set the energy radiated away by the particles
    void setNrjRadiation( double value )
    {
        nrj_radiation = value;
    }

    //! Add the energy radiated away by the particles
    void addNrjRadiation( double value )
    {
        nrj_radiation += value;
    }

    //! Set gained via new particles
    void setNewParticlesNRJ( double value )
    {
        new_particles_energy_ = value;
    }
    
    //! Get energy gained via new particles
    double getNewParticlesNRJ() const
    {
        return new_particles_energy_;
    }

    //! Reinitialize the scalar diagnostics buffer
    void reinitDiags()
    {
        //nrj_bc_lost = 0;
        nrj_mw_lost = 0;
        //new_particles_energy_ = 0;
        //nrj_radiation = 0;
    }

    inline void storeNRJlost( double nrj )
    {
        nrj_mw_lost += nrj;
    };

    inline double computeNRJ()
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

    inline int getMemFootPrint()
    {
        /*int speciesSize  = ( 2*nDim_particle + 3 + 1 )*sizeof(double) + sizeof(short);
        if ( particles->is_test )
            speciesSize += sizeof ( unsigned int );*/
        //speciesSize *= getNbrOfParticles();
        int speciesSize( 0 );
        speciesSize += particles->double_prop.size()*sizeof( double );
        speciesSize += particles->short_prop.size()*sizeof( short );
        speciesSize += particles->uint64_prop.size()*sizeof( uint64_t );
        speciesSize *= getParticlesCapacity();
        return speciesSize;
    }

    //! Method to import particles in this species while conserving the sorting among bins
    virtual void importParticles( Params &, Patch *, Particles &, std::vector<Diagnostic *> & );

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

protected:

    //! Patch length
    unsigned int length_[3];

private:
    //! Number of steps for Maxwell-Juettner cumulative function integration
    //! \todo{Put in a code constant class}

//    unsigned int nE;

    //! Parameter used when defining the Maxwell-Juettner function (corresponds to a Maximum energy)
//    double muEmax;

    //! Parameter used when defining the Maxwell-Juettner function (corresponds to a discretization step in energy)
//    double dE;

};

#endif
