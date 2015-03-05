#ifndef PATCH_H
#define PATCH_H

#include <vector>

#include "SpeciesFactory.h"
#include "ElectroMagnFactory.h"
#include "InterpolatorFactory.h"
#include "ProjectorFactory.h"

#include "PicParams.h"
#include "LaserParams.h"
#include "SmileiMPI.h"


//! Class Patch : sub MPI domain 
//!     Collection of patch = MPI domain
class Patch
{

public:
    //! Constructor for Patch
    Patch(PicParams& params, LaserParams& laser_params, SmileiMPI* smpi) {

	vecSpecies = SpeciesFactory::createVector(params, smpi);              // + patchId + min_loc/cell_index(ref smpi, creta Pos & sort) + new n_space
	// -> partBoundCond : min/max_loc (smpi)
	EMfields   = ElectroMagnFactory::create(params, laser_params, smpi);  // + patchId + new n_space (now = params by smpi) + BC
	Interp     = InterpolatorFactory::create(params, smpi);               // + patchId -> idx_domain_begin (now = ref smpi)
	Proj       = ProjectorFactory::create(params, smpi);                  // + patchId -> idx_domain_begin (now = ref smpi)
	
    };

    //! Destructor for Patch
    ~Patch() {

	delete Proj;
	delete Interp;
	delete EMfields;
	for (unsigned int ispec=0 ; ispec<vecSpecies.size(); ispec++) delete vecSpecies[ispec];
	vecSpecies.clear();
	    
    };

    std::vector<Species*> vecSpecies;
    ElectroMagn* EMfields;

    Interpolator* Interp;
    Projector* Proj;

protected:

private:
    //! cell_starting_global_index : index of 1st cell of local sub-subdomain in the global domain
    //!     - concerns ghost data
    //!     - "- oversize" on rank 0
    std::vector<int> cell_starting_global_index;
    //! "Real" min limit of local sub-subdomain (ghost data not concerned)
    //!     - "0." on rank 0
    std::vector<double> min_local;
    //! "Real" max limit of local sub-subdomain (ghost data not concerned)
    std::vector<double> max_local;

    //! number of cells in every direction of the local sub-subdomain, all patch have the same size ???
    std::vector<unsigned int> n_space;

};

#endif
