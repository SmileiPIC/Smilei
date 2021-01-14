#ifndef SPECIESTOMP_H
#define SPECIESTOMP_H

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
class Species_taskomp: public Species
{
public:
    // -----------------------------------------------------------------------------
    //  Constructor/destructor

    //! Species constructor
    Species_taskomp( Params &, Patch * );

    //! Species destructor
    virtual ~Species_taskomp();

    
    //! Method calculating the Particle dynamics (interpolation, pusher, projection)
    void dynamics( double time, unsigned int ispec,
                           ElectroMagn *EMfields,
                           Params &params, bool diag_flag,
                           PartWalls *partWalls, Patch *patch, SmileiMPI *smpi,
                           RadiationTables &RadiationTables,
                           MultiphotonBreitWheelerTables &MultiphotonBreitWheelerTables,
                           std::vector<Diagnostic *> &localDiags ) override;
    //! Method calculating the Particle dynamics (interpolation, pusher, projection)
    void dynamicsWithTasks( double time, unsigned int ispec,
                           ElectroMagn *EMfields,
                           Params &params, bool diag_flag,
                           PartWalls *partWalls, Patch *patch, SmileiMPI *smpi,
                           RadiationTables &RadiationTables,
                           MultiphotonBreitWheelerTables &MultiphotonBreitWheelerTables,
                           std::vector<Diagnostic *> &localDiags, int buffer_id );

    
    int *bin_has_pushed;
    int *bin_has_done_particles_BC;
    int *bin_has_projected;

    std::vector<double *> b_Jx;
    std::vector<double *> b_Jy;
    std::vector<double *> b_Jz;
    std::vector<double *> b_rho;

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
