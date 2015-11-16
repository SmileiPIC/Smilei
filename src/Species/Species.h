#ifndef SPECIES_H
#define SPECIES_H

#include <vector>
#include <string>

#include "Particles.h"
#include "Params.h"
#include "Pusher.h"
//#include "PartBoundCond.h"

#include "Params.h"
#include "Pusher.h"
#include "Ionization.h"
#include "ElectroMagn.h"
#include "Profile.h"
#include "SpeciesMPI.h"

class ElectroMagn;
class Pusher;
class Interpolator;
class Projector;
class PartBoundCond;
class PartWall;
class Field3D;
class Patch;

//! class Species
class Species
{
public:
    SpeciesMPI specMPI;

    //! Species creator
    Species(Params&, int, Patch*);

    void initCluster(Params&);
    void initSpecies(Params&);

    //! Species destructor
    virtual ~Species();

    //! Species index
    unsigned int speciesNumber;

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
    // TO do defince particles.capacity = min.capacity
    inline unsigned int getParticlesCapacity() const {
        return (*particles).capacity();
    }

    //! Method calculating the Particle dynamics (interpolation, pusher, projection)
    virtual void dynamics(double time, unsigned int ispec, ElectroMagn* EMfields, Interpolator* interp,
                          Projector* proj, Params &params, int diag_flag,
			  std::vector<PartWall*> vecPartWall, Patch* patch);

    //! Method used to initialize the Particle position in a given cell
    void initPosition(unsigned int, unsigned int, double *, unsigned int, std::vector<double>, std::string);

    //! Method used to initialize the Particle 3d momentum in a given cell
    void initMomentum(unsigned int, unsigned int, double *, double *, std::string, std::vector<double>&);

    //! Method used to initialize the Particle weight (equivalent to a charge density) in a given cell
    void initWeight(unsigned int, unsigned int, unsigned int, double);

    //! Method used to initialize the Particle charge
    void initCharge(unsigned int, unsigned int, unsigned int, double);
    
    //! Maximum charge at initialization
    double max_charge;

    //! Method used to sort particles
    void sort_part();
    void count_sort_part(Params& param);

    void updateMvWinLimits(double x_moved);

    //! Vector containing all Particles of the considered Species
    Particles *particles;
    Particles particles_sorted[2];
    //std::vector<int> index_of_particles_to_exchange;

    //! Ionization method
    Ionization* Ionize;

    //! to keep rack of ionized electrons
    Species *electron_species;

    //! Cluster width in number of cells
    unsigned int clrw; //Should divide the number of cells in X of a single MPI domain. 
    //! first and last index of each particle bin
    std::vector<int> bmin, bmax;
    //! sub dimensions of buffers for dim > 1
    unsigned int b_dim0, b_dim1, b_dim2, b_lastdim;

    //! Size of the projection buffer
    unsigned int size_proj_buffer;

    //! Oversize (copy from Params)
    std::vector<unsigned int> oversize;

    //! MPI structure to exchange particles
    MPI_Datatype typePartSend ;

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

    //Copy of the species parameters from Params
    SpeciesStructure species_param;

    //! Method to know if we have to project this species or not.
    bool  isProj(double time_dual, SimWindow* simWindow);

    double getLostNrjBC() const {return species_param.mass*nrj_bc_lost;}
    double getLostNrjMW() const {return species_param.mass*nrj_mw_lost;}

    double getNewParticlesNRJ() const {return species_param.mass*nrj_new_particles;}
    void reinitDiags() { 
	nrj_bc_lost = 0;
	nrj_mw_lost = 0;
	nrj_new_particles = 0;
    }
    inline void storeNRJlost( double nrj ) { nrj_mw_lost = nrj; };

    inline double computeNRJ() {
	double nrj(0.);
	for ( unsigned int iPart=0 ; iPart<getNbrOfParticles() ; iPart++ )
	    nrj += (*particles).weight(iPart)*((*particles).lor_fac(iPart)-1.0);
	return nrj;
    }

    inline int getMemFootPrint() {
	int speciesSize  = ( 2*ndim + 3 + 1 )*sizeof(double) + sizeof(short);
	if ( particles->isTestParticles )
	    speciesSize += sizeof ( unsigned int );
	//speciesSize *= getNbrOfParticles();
	speciesSize *= getParticlesCapacity();
	return speciesSize;
    }

private:
    
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
    
    Profile *ppcProfile;
    
    //! 2 times pi
    double PI2;
    double dx_inv_, dy_inv_;
    
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


    //! sub primal dimensions of fields
    unsigned int f_dim0, f_dim1, f_dim2;

    //! Method used to apply boundary-condition for the Particles of the considered Species
    PartBoundCond* partBoundCond;
    
    //! Method used to Push the particles (change momentum & change position)
    Pusher* Push;

    //! Method to create new particles.
    int  createParticles(std::vector<unsigned int> n_space_to_create, std::vector<double> cell_index, int new_bin_idx,  Params& param);

    //! Accumulate nrj lost with bc
    double nrj_bc_lost;
    //! Accumulate nrj lost with moving window
    double nrj_mw_lost;
    //! Accumulate nrj added with new particles
    double nrj_new_particles;
};

#endif
