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
	unsigned int npatches, firstpatch;
        unsigned int m0, m1, m2; //Defines the log2 of the total number of patches for each direction.
        m0 = 3;
        m1 = 2;
        m2 = 0;
        npatches = (1 << (m0 + m1 + m2)) / smpi->getSize() + ( smpi->getRank() < (1 << (m0 + m1 + m2))%smpi->getSize() );// npatches = 2^(m0+m1+m2) / number of mpi process. Local number of patches. 
        //Naive initialization of patch_count, assuming all mpi processes initially have the same number of patches.
        for (unsigned int impi = 0 ; impi < smpi->getSize() ; impi++) {
            smpi->patch_count[impi] = (1 << (m0 + m1 + m2)) / smpi->getSize() +1*( impi < (1 << (m0 + m1 + m2))%smpi->getSize() );
        }
        firstpatch = 0;
        for (unsigned int impi = 0 ; impi < smpi->getRank() ; impi++) {
            firstpatch += smpi->patch_count[impi];
        }
        

	// Modified to test Patch integration
	//npatches = m0*m1;
        vecPatches.resize(npatches);
	//std::cout << "n_space : " << params.n_space[0] << " " << params.n_space[1] << std::endl;
	//params.n_space[0] /= m0;
	//params.n_space[1] /= m1;
	//std::cout << "n_space : " << params.n_space[0] << " " << params.n_space[1] << std::endl;


        // create species
        unsigned int nPart;
        for (unsigned int ipatch = 0 ; ipatch < npatches ; ipatch++) {
	    // For Patch integration : ipatch is local to MPI process
	    vecPatches[ipatch] = PatchesFactory::create(params, laser_params, smpi, m0, m1, m2, firstpatch + ipatch);
            //vecPatches[ipatch] = PatchesFactory::create(params, laser_params, smpi, m0, m1, m2, ipatch);
        }

        return vecPatches;
    }

};

#endif
