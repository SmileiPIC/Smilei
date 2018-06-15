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
    SpeciesV(Params&, Patch*);
    //! Species destructor
    virtual ~SpeciesV();

    void initCluster(Params& params) override;

//virtual void dynamics(double time, unsigned int ispec, ElectroMagn* EMfields, Interpolator* interp,
//                      Projector* proj, Params &params, bool diag_flag,
//                      PartWalls* partWalls, Patch* patch, SmileiMPI* smpi) override;
//

    //! Method calculating the Particle dynamics (interpolation, pusher, projection)
    void dynamics(double time, unsigned int ispec,
                          ElectroMagn* EMfields,
                          Interpolator* interp,
                          Projector* proj, Params &params, bool diag_flag,
                          PartWalls* partWalls, Patch* patch, SmileiMPI* smpi,
                          RadiationTables &RadiationTables,
                          MultiphotonBreitWheelerTables & MultiphotonBreitWheelerTables,
                          std::vector<Diagnostic*>& localDiags) override;

    //! Method calculating the Particle updated momentum (interpolation, momentum pusher, only particles interacting with envelope)
    void ponderomotive_update_susceptibility_and_momentum(double time_dual, unsigned int ispec,
                               ElectroMagn* EMfields, Interpolator* Interp_envelope, Projector* Proj_susceptibility,
                               Params &params, bool diag_flag,
                               Patch* patch, SmileiMPI* smpi,
                               std::vector<Diagnostic*>& localDiags) override;

    //! Method calculating the Particle updated position (interpolation, position pusher, only particles interacting with envelope)
    // and projecting charge density and thus current density (through Esirkepov method) for Maxwell's Equations
    void ponderomotive_update_position_and_currents(double time_dual, unsigned int ispec,
                               ElectroMagn* EMfields, Interpolator* Interp_envelope, Projector* Proj,
                               Params &params, bool diag_flag, PartWalls* partWalls,
                               Patch* patch, SmileiMPI* smpi,
                               std::vector<Diagnostic*>& localDiags) override;

    //! Method calculating the Particle charge on the grid (projection)
    void computeCharge(unsigned int ispec, ElectroMagn* EMfields, Projector* Proj) override;

    //! Method used to sort particles
    void sort_part(Params& params) override;
    //void count_sort_part(Params& param);
    void compute_part_cell_keys(Params &params);
    
    //! 
    void add_space_for_a_particle() override {
        particles->cell_keys.push_back(-1);
    }

    //! Method to import particles in this species while conserving the sorting among bins
    void importParticles( Params&, Patch*, Particles&, std::vector<Diagnostic*>& )override;

private:

};

#endif

