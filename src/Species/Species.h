#ifndef SPECIES_H
#define SPECIES_H

#include <vector>
#include <string>

#include "Particles.h"
#include "PicParams.h"
#include "Pusher.h"
//#include "PartBoundCond.h"
#include "PicParams.h"
#include "SmileiMPI.h"
#include "Pusher.h"
#include "Ionization.h"
#include "ElectroMagn.h"
#include "DensityProfile.h"

class ElectroMagn;
class Pusher;
class Interpolator;
class Projector;
class PartBoundCond;
class Field3D;

//! class Species
class Species
{
public:
    //! Species creator
    Species(PicParams&, int, SmileiMPI*);

    //! Species destructor
    virtual ~Species();

    //! Species index
    unsigned int speciesNumber;

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
    // TO do defince particles.capacity = min.capacity
    inline unsigned int getParticlesCapacity() const {
        return particles.capacity();
    }

    //! Method calculating the Particle dynamics (interpolation, pusher, projection)
    virtual void dynamics(double time, unsigned int ispec, ElectroMagn* EMfields, Interpolator* interp,
                          Projector* proj, SmileiMPI *smpi, PicParams &params, SimWindow* simWindow);

    //! Method used to initialize the Particle position in a given cell
    void initPosition(unsigned int, unsigned int, double *, unsigned int, std::vector<double>, std::string);

    //! Method used to initialize the Particle 3d momentum in a given cell
    void initMomentum(unsigned int, unsigned int, double *, double *, std::string, std::vector<double>&);

    //! Method used to initialize the Particle weight (equivalent to a charge density) in a given cell
    void initWeight(PicParams*, unsigned int, unsigned int, double);

    //! Method used to initialize the Particle charge
    void initCharge(PicParams*, unsigned int, unsigned int, double);

    //! Method used to save all Particles properties for the considered Species
    void dump(std::ofstream&);

    //! Method used to sort particles
    void sort_part();

    void movingWindow_x(unsigned int shift, SmileiMPI *smpi, PicParams& param);
    void defineNewCells(unsigned int shift, SmileiMPI *smpi, PicParams& param);
    void updateMvWinLimits(double x_moved);

    //! Vector containing all Particles of the considered Species
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

    //! Oversize (copy from picparams)
    std::vector<unsigned int> oversize;

    //! Cell_length (copy from picparams)
    std::vector<double> cell_length;

    inline void clearExchList(int tid) {
	    indexes_of_particles_to_exchange_per_thd[tid].clear();
    }
    inline void addPartInExchList(int tid, int iPart) {
        indexes_of_particles_to_exchange_per_thd[tid].push_back(iPart);
    }
    std::vector< std::vector<int> > indexes_of_particles_to_exchange_per_thd;

    //Copy of the species parameters from picparams
    SpeciesStructure species_param;

    //! Method to know if we have to project this species or not.
    bool  isProj(double time_dual, SimWindow* simWindow);

    double getLostNrjBC() const {return species_param.mass*nrj_bc_lost;}
    double getLostNrjMW() const {return species_param.mass*nrj_mw_lost;}

    double getNewParticlesNRJ() const {return species_param.mass*nrj_new_particles;}

private:
    
    //! vector of density (one per species)
    DensityProfile *densityProfile;
    
    //! 2 times pi
    double PI2;
    
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

    //! Method used to apply boundary-condition for the Particles of the considered Species
    PartBoundCond* partBoundCond;

    //! Method used to Push the particles (change momentum & change position)
    Pusher* Push;

    //! Method to create new particles.
    int  createParticles(std::vector<unsigned int> n_space_to_create, std::vector<double> cell_index, int new_bin_idx,  PicParams& param);

    //! Accumulate nrj lost with bc
    double nrj_bc_lost;
    //! Accumulate nrj lost with moving window
    double nrj_mw_lost;
    //! Accumulate nrj added with new particles during last call to createParticles
    double nrj_new_particles;
};

#endif
