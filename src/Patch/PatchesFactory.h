#ifndef PATCHESFACTORY_H
#define PATCHESFACTORY_H

#include "Patch.h"

#include "Tools.h"

class PatchesFactory {
public:
    static Patch* create(PicParams& params, LaserParams& laser_params, SmileiMPI* smpi) {
        Patch* patch = new Patch(params, laser_params, smpi);
        return patch;
    }

    static std::vector<Patch*> createVector(PicParams& params, LaserParams& laser_params, SmileiMPI* smpi) {
        std::vector<Patch*> vecPatches;

	// Compute npatches (1 is std MPI behavior)
	int npatches(1);
        vecPatches.resize(npatches);

        // create species
        unsigned int nPart;
        for (unsigned int ipatch=0 ; ipatch<npatches ; ipatch++) {
            vecPatches[ipatch] = PatchesFactory::create(params, laser_params, smpi);
        }

        return vecPatches;
    }

};

#endif
