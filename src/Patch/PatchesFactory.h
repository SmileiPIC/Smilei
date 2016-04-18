#ifndef PATCHESFACTORY_H
#define PATCHESFACTORY_H

#include "VectorPatch.h"
#include "Patch1D.h"
#include "Patch2D.h"

#include "DiagsVectorPatch.h"

#include "Tools.h"

class PatchesFactory {
public:
    
    // Create one patch from scratch
    static Patch* create(Params& params, SmileiMPI* smpi, unsigned int ipatch, unsigned int n_moved=0) {
        Patch* patch;
        if (params.geometry == "1d3v")
            patch = new Patch1D(params, smpi, ipatch, n_moved);
        else 
            patch = new Patch2D(params, smpi, ipatch, n_moved);
        return patch;
    }
    
    // Clone one patch (avoid reading again the namelist)
    static Patch* clone(Patch* patch, Params& params, SmileiMPI* smpi, unsigned int ipatch, unsigned int n_moved=0) {
        Patch* newPatch;
        if (params.geometry == "1d3v")
            newPatch = new Patch1D(params, smpi, ipatch, n_moved);
        else 
            newPatch = new Patch2D(params, smpi, ipatch, n_moved);
        return newPatch;
    }
    
    // Create a vector of patches
    static VectorPatch createVector(Params& params, SmileiMPI* smpi) {
        VectorPatch vecPatches;
        
        // Compute npatches (1 is std MPI behavior)
        unsigned int npatches, firstpatch;
        npatches = smpi->patch_count[smpi->getRank()];// Number of patches owned by current MPI process.
        firstpatch = 0;
        for (unsigned int impi = 0 ; impi < smpi->getRank() ; impi++) {
            firstpatch += smpi->patch_count[impi];
        }
#ifdef _DEBUGPATCH
        std::cout << smpi->getRank() << ", nPatch = " << npatches << " - starting at " << firstpatch << std::endl;        
#endif
        
        // create patches (create patch#0 then clone it)
        vecPatches.resize(npatches);
        vecPatches.patches_[0] = create(params, smpi, firstpatch);
        for (unsigned int ipatch = 1 ; ipatch < npatches ; ipatch++) {
            vecPatches.patches_[ipatch] = clone(vecPatches(0), params, smpi, firstpatch + ipatch);
        }
        vecPatches.set_refHindex();
        
        // Patch initializations which needs some sync (parallel output, are data distribution)
        int itime(0);
        DiagsVectorPatch::initDumpFields(vecPatches, params, itime);
        DiagsVectorPatch::initCollisions(vecPatches, params, smpi );
        
        vecPatches.update_field_list();
        
        // Figure out if there are antennas
        vecPatches.hasAntennas = ( vecPatches(0)->EMfields->antennas.size() > 0 );
        
        vecPatches.createGlobalDiags( params, smpi );
        vecPatches.initAllDiags( params, smpi );
        
        return vecPatches;
    }

};

#endif
