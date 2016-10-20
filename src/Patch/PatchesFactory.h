#ifndef PATCHESFACTORY_H
#define PATCHESFACTORY_H

#include "VectorPatch.h"
#include "Patch1D.h"
#include "Patch2D.h"
#include "Patch3D.h"

#include "Tools.h"

class PatchesFactory {
public:
    
    // Create one patch from scratch
    static Patch* create(Params& params, SmileiMPI* smpi, unsigned int ipatch, unsigned int n_moved=0) {
        Patch* patch;
        if (params.geometry == "1d3v")
            patch = new Patch1D(params, smpi, ipatch, n_moved);
        else if (params.geometry == "2d3v") 
            patch = new Patch2D(params, smpi, ipatch, n_moved);
        else if (params.geometry == "3d3v") 
            patch = new Patch3D(params, smpi, ipatch, n_moved);
        return patch;
    }
    
    // Clone one patch (avoid reading again the namelist)
    static Patch* clone(Patch* patch, Params& params, SmileiMPI* smpi, unsigned int ipatch, unsigned int n_moved=0, bool with_particles = true) {
        Patch* newPatch;
        if (params.geometry == "1d3v")
            newPatch = new Patch1D(static_cast<Patch1D*>(patch), params, smpi, ipatch, n_moved, with_particles);
        else if (params.geometry == "2d3v")
            newPatch = new Patch2D(static_cast<Patch2D*>(patch), params, smpi, ipatch, n_moved, with_particles);
        else if (params.geometry == "3d3v")
            newPatch = new Patch3D(static_cast<Patch3D*>(patch), params, smpi, ipatch, n_moved, with_particles);
        return newPatch;
    }
    
    // Create a vector of patches
    static VectorPatch createVector(Params& params, SmileiMPI* smpi) {
        VectorPatch vecPatches;
        
        // Compute npatches (1 is std MPI behavior)
        unsigned int npatches, firstpatch;
        npatches = smpi->patch_count[smpi->getRank()];// Number of patches owned by current MPI process.
        firstpatch = 0;
        for (unsigned int impi = 0 ; impi < (unsigned int)smpi->getRank() ; impi++) {
            firstpatch += smpi->patch_count[impi];
        }
        DEBUG( smpi->getRank() << ", nPatch = " << npatches << " - starting at " << firstpatch );
        
        // Create patches (create patch#0 then clone it)
        vecPatches.resize(npatches);
        vecPatches.patches_[0] = create(params, smpi, firstpatch);
        
        TITLE("Initializing Patches");
        MESSAGE(1,"First patch created");
        unsigned int percent=10;
        for (unsigned int ipatch = 1 ; ipatch < npatches ; ipatch++) {
            if( (100*ipatch)/npatches > percent ) {
                MESSAGE(2,"Approximately "<<percent<<"% of patches created");
                percent += 10;
            }
            vecPatches.patches_[ipatch] = clone(vecPatches(0), params, smpi, firstpatch + ipatch);
        }
        MESSAGE(1,"All patches created");
        // print number of particles
        vecPatches.printNumberOfParticles( smpi );
        
        vecPatches.set_refHindex();
        
        vecPatches.update_field_list();
        
        TITLE("Initializing Diagnostics, antennas, and external fields")
        vecPatches.createDiags( params, smpi );
        
        // Figure out if there are antennas
        vecPatches.nAntennas = vecPatches(0)->EMfields->antennas.size();
        vecPatches.initExternals( params );
        
        MESSAGE(1,"Done initializing diagnostics");
        return vecPatches;
    }

};

#endif
