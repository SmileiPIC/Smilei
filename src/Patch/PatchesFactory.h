#ifndef PATCHESFACTORY_H
#define PATCHESFACTORY_H

#include "Patch.h"

#include "Tools.h"

class PatchesFactory {
public:
    static Patch* create(PicParams& params, LaserParams& laser_params, SmileiMPI* smpi, unsigned int m0,unsigned int  m1,unsigned int  m2,unsigned int  ipatch) {
        Patch* patch = new Patch(params, laser_params, smpi, m0, m1, m2, ipatch);
        return patch;
    }

    static std::vector<Patch*> createVector(PicParams& params, LaserParams& laser_params, SmileiMPI* smpi) {
        std::vector<Patch*> vecPatches;

	// Compute npatches (1 is std MPI behavior)
	unsigned int npatches;
        unsigned int m0, m1, m2; //Defines the total number of patches.
        m0 = 0;
        m1 = 0;
        m2 = 0;
        npatches = (1 << (m0 + m1 + m2)) / smpi->getSize() ;// npatches = 2^(m0+m1+m2) / number of mpi process. Local number of patches. 
        vecPatches.resize(npatches);


        // create species
        unsigned int nPart;
        for (unsigned int ipatch = smpi->getRank()*npatches ; ipatch<(smpi->getRank()+1)*npatches ; ipatch++) {
            vecPatches[ipatch] = PatchesFactory::create(params, laser_params, smpi, m0, m1, m2, ipatch);
        }

        return vecPatches;
    }

};

#endif
