#ifndef PATCHESFACTORY_H
#define PATCHESFACTORY_H

#include "Patch.h"

#include "Tools.h"

class PatchesFactory {
public:
    static Patch* create(PicParams& params, DiagParams& diag_params, LaserParams& laser_params, SmileiMPI* smpi, unsigned int  ipatch) {
	Patch* patch = new Patch(params, diag_params, laser_params, smpi, ipatch);
        return patch;
    }

    static VectorPatch createVector(PicParams& params, DiagParams& diag_params, LaserParams& laser_params, SmileiMPI* smpi) {
        VectorPatch vecPatches;

	// Compute npatches (1 is std MPI behavior)
	unsigned int npatches, firstpatch;
        //unsigned int m0, m1, m2; //Defines the log2 of the total number of patches for each direction.
        //m0 = 0;
        //m1 = 0;
        //m2 = 0;
	////std::cout << " params.number_of_patches[0] = " << params.number_of_patches[0] << std::endl;
        //while ((params.number_of_patches[0] >> m0) >1) m0++ ;
	////std::cout << " m0 = " << m0 << std::endl;
        //while ((params.number_of_patches[1] >> m1) >1) m1++ ;
        //while ((params.number_of_patches[2] >> m2) >1) m2++ ;
        npatches = smpi->patch_count[smpi->getRank()];// Number of patches owned by current MPI process.
        firstpatch = 0;
        for (unsigned int impi = 0 ; impi < smpi->getRank() ; impi++) {
            firstpatch += smpi->patch_count[impi];
        }
	std::cout << smpi->getRank() << ", nPatch = " << npatches << " - starting at " << firstpatch << std::endl;        

	// Modified to test Patch integration

	std::cout << "n_space : " << params.n_space[0] << " " << params.n_space[1] << std::endl;
	std::cout << "n_patch : " << params.number_of_patches[0] << " " << params.number_of_patches[1] << std::endl;

	params.n_space[0] /= params.number_of_patches[0];
	params.n_space[1] /= params.number_of_patches[1];
	std::cout << "\ Patch : n_space : " << params.n_space[0] << " " << params.n_space[1] << std::endl;


        // create species
        vecPatches.resize(npatches);
        for (unsigned int ipatch = 0 ; ipatch < npatches ; ipatch++) {
	    vecPatches.patches_[ipatch] = PatchesFactory::create(params, diag_params, laser_params, smpi, firstpatch + ipatch);
        }

        return vecPatches;
    }

};

#endif
