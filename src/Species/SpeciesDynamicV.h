#ifndef SPECIESDYNAMICV_H
#define SPECIESDYNAMICV_H

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
class SpeciesDynamicV : public SpeciesV
{
 public:
    //! Species creator
    SpeciesDynamicV(Params&, Patch*);
    //! Species destructor
    virtual ~SpeciesDynamicV();

    void resizeCluster(Params& params) override;

    //! This function configures the species according to the vectorization mode
    void configuration( Params& params, Patch * patch) override;

    //! This function reconfigures the species according to the vectorization mode
    void reconfiguration( Params& params, Patch * patch) override;

    //! This function reconfigures the species operators
    void reconfigure_operators(Params& param, Patch  * patch);

    //! This function reconfigures the species to be imported
    void reconfigure_particle_importation();

    //! Compute cell_keys for all particles of the current species
    void compute_part_cell_keys(Params &params);

    //! Method to import particles in this species while conserving the sorting among bins
    void importParticles( Params&, Patch*, Particles&, std::vector<Diagnostic*>& )override;

private:

    // Metrics for the dynamic vectorization
    int max_number_of_particles_per_cells;
    int min_number_of_particles_per_cells;
    double ratio_number_of_vecto_cells;

    //! Number of packs of particles that divides the total number of particles
    unsigned int npack_;
    //! Size of the pack in number of particles
    unsigned int packsize_;

};

#endif
