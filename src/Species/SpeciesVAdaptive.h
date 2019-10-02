#ifndef SPECIESADAPTIVE_H
#define SPECIESADAPTIVE_H

#include <vector>
#include <string>

#include "SpeciesV.h"

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
class SpeciesVAdaptive : public SpeciesV
{
public:
    //! Species creator
    SpeciesVAdaptive( Params &, Patch * );
    //! Species destructor
    virtual ~SpeciesVAdaptive();
    
    //! Method calculating the Particle dynamics (interpolation, pusher, projection)
    //! without vectorized operators but with the cell sorting algorithm
    void scalarDynamics( double time, unsigned int ispec,
                          ElectroMagn *EMfields,
                          Params &params, bool diag_flag,
                          PartWalls *partWalls, Patch *patch, SmileiMPI *smpi,
                          RadiationTables &RadiationTables,
                          MultiphotonBreitWheelerTables &MultiphotonBreitWheelerTables,
                          std::vector<Diagnostic *> &localDiags ) override;
                          
    //! This function configures the type of species according to the default mode
    //! regardless the number of particles per cell
    void defaultConfigure( Params &params, Patch *patch ) override;
    
    //! This function configures the species according to the vectorization mode
    void configuration( Params &params, Patch *patch ) override;
    
    //! This function reconfigures the species according to the vectorization mode
    void reconfiguration( Params &params, Patch *patch ) override;
    
    //! This function reconfigures the species operators
    void reconfigure_operators( Params &param, Patch   *patch );
    
    //void count_sortParticles(Params& param);
    //void computeParticleCellKeys(Params &params);
    
    void scalarPonderomotiveUpdateSusceptibilityAndMomentum( double time_dual, unsigned int ispec,
            ElectroMagn *EMfields,
            Params &params, bool diag_flag,
            Patch *patch, SmileiMPI *smpi,
            std::vector<Diagnostic *> &localDiags ) override;
            
    void scalarPonderomotiveUpdatePositionAndCurrents( double time_dual, unsigned int ispec,
            ElectroMagn *EMfields,
            Params &params, bool diag_flag, PartWalls *partWalls,
            Patch *patch, SmileiMPI *smpi,
            std::vector<Diagnostic *> &localDiags ) override;
            
            
private:

    // Metrics for the adaptive vectorization
    //int max_number_of_particles_per_cells;
    //int min_number_of_particles_per_cells;
    //double ratio_number_of_vecto_cells;
    
    //! Number of packs of particles that divides the total number of particles
    unsigned int npack_;
    //! Size of the pack in number of particles
    unsigned int packsize_;
    
};

#endif
