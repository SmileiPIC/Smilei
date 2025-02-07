#ifndef SPECIESV_H
#define SPECIESV_H

#include <vector>
#include <string>

#include "Species.h"

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
class SpeciesV : public Species
{
public:
    //! Species creator
    SpeciesV( Params &, Patch * );
    //! Species destructor
    virtual ~SpeciesV();

    void initCluster( Params &params, Patch *patch ) override;
    
    //! Method calculating the Particle dynamics (interpolation, pusher, projection)
    void dynamics( double time, unsigned int ispec,
                   ElectroMagn *EMfields,
                   Params &params, bool diag_flag,
                   PartWalls *partWalls, Patch *patch, SmileiMPI *smpi,
                   RadiationTables &RadiationTables,
                   MultiphotonBreitWheelerTables &MultiphotonBreitWheelerTables ) override;

    //! Method projecting susceptibility and calculating the particles updated momentum (interpolation, momentum pusher), only particles interacting with envelope
    void ponderomotiveUpdateSusceptibilityAndMomentum( double time_dual, 
            ElectroMagn *EMfields,
            Params &params, 
            Patch *patch, SmileiMPI *smpi ) override;

    //! Method projecting susceptibility, only particles interacting with envelope
    void ponderomotiveProjectSusceptibility( double time_dual,
            ElectroMagn *EMfields,
            Params &params, 
            Patch *patch, SmileiMPI *smpi ) override;


    //! Method calculating the Particle updated position (interpolation, position pusher, only particles interacting with envelope)
    // and projecting charge density and thus current density (through Esirkepov method) for Maxwell's Equations
    void ponderomotiveUpdatePositionAndCurrents( double time_dual, unsigned int ispec,
            ElectroMagn *EMfields,
            Params &params, bool diag_flag, PartWalls *partWalls,
            Patch *patch, SmileiMPI *smpi ) override;

    //! Method calculating the Particle charge on the grid (projection)
    void computeCharge( ElectroMagn *EMfields, bool old=false ) override;

    //! Method used to sort particles
    void sortParticles( Params &params ) override;
    //void countSortParticles(Params& param);

    //! Compute cell_keys for all particles from istart to iend
    void computeParticleCellKeys(   Params    & params,
                                    Particles * particles,
                                    int       * __restrict__ cell_keys,
                                    int       * __restrict__ count,
                                    unsigned int istart,
                                    unsigned int iend ) override;

    //! Compute cell_keys for all particles of the current species
    void computeParticleCellKeys( Params &params ) override;

    //! Create a new entry for a particle
    void addSpaceForOneParticle() override
    {
        particles->cell_keys.push_back( -1 );
    }

    //! Method to import particles in this species while conserving the sorting among bins
    void importParticles( Params &, Patch *, Particles &, std::vector<Diagnostic *> &, double, Ionization *I = nullptr )override;

    //! Method performing the merging of particles
    virtual void mergeParticles( double time_dual )override;


private:

    //! Number of packs of particles that divides the total number of particles
    unsigned int npack_;
    //! Size of the pack in number of particles
    unsigned int packsize_;

};

#endif
