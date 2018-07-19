#ifndef SPECIESDYNAMICV2_H
#define SPECIESDYNAMICV2_H

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
class SpeciesDynamicV2 : public SpeciesV
{
 public:
    //! Species creator
    SpeciesDynamicV2(Params&, Patch*);
    //! Species destructor
    virtual ~SpeciesDynamicV2();

    //! Method calculating the Particle dynamics (interpolation, pusher, projection)
    //! without vectorized operators but with the cell sorting algorithm
    void scalar_dynamics(double time, unsigned int ispec,
                        ElectroMagn* EMfields,
                        Params &params, bool diag_flag,
                        PartWalls* partWalls, Patch* patch, SmileiMPI* smpi,
                        RadiationTables &RadiationTables,
                        MultiphotonBreitWheelerTables & MultiphotonBreitWheelerTables,
                        std::vector<Diagnostic*>& localDiags);

    //! This function configures the species according to the vectorization mode
    void configuration( Params& params, Patch * patch) override;

    //! This function reconfigures the species according to the vectorization mode
    void reconfiguration( Params& params, Patch * patch) override;

    //! This function reconfigures the species operators
    void reconfigure_operators(Params& param, Patch  * patch);

    //void count_sort_part(Params& param);
    //void compute_part_cell_keys(Params &params);

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
