#ifndef PATCHESFACTORY_H
#define PATCHESFACTORY_H

#include "VectorPatch.h"
#include "Patch1D.h"
#include "Patch2D.h"

#include "DiagsVectorPatch.h"
#include "SpeciesFactory.h"
#include "ElectroMagnFactory.h"
#include "InterpolatorFactory.h"
#include "ProjectorFactory.h"
#include "SmileiIOFactory.h"
#include "DiagnosticFactory.h"

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
        
        // First, do a series of initializations
        patch->initStep1(params);
        patch->initStep2(params); // depends on 1D, 2D, ...
        patch->initStep3(params, smpi, n_moved);
        
        // initialize vector of Species (virtual)
        patch->vecSpecies = SpeciesFactory::createVector(params, patch);
        
        // initialize the electromagnetic fields (virtual)
        patch->EMfields   = ElectroMagnFactory::create(params, patch->vecSpecies, patch);
        
        // interpolation operator (virtual)
        patch->Interp     = InterpolatorFactory::create(params, patch);               // + patchId -> idx_domain_begin (now = ref smpi)
        // projection operator (virtual)
        patch->Proj       = ProjectorFactory::create(params, patch);                  // + patchId -> idx_domain_begin (now = ref smpi)
        
        // Initialize local diags
        patch->localDiags = DiagnosticFactory::createLocalDiagnostics(params, smpi, patch);
        patch->sio = SmileiIOFactory::create(params, patch);
        
        // Initialize the collisions
        patch->vecCollisions = Collisions::create(params, patch, patch->vecSpecies);
        
        // Initialize the particle walls
        patch->partWalls = new PartWalls(params, patch);

        patch->createType(params);
        return patch;
    }
    
    // Clone one patch (avoid reading again the namelist)
    static Patch* clone(Patch* patch, Params& params, SmileiMPI* smpi, unsigned int ipatch, unsigned int n_moved=0) {
        Patch* newPatch;
        if (params.geometry == "1d3v")
            newPatch = new Patch1D(params, smpi, ipatch, n_moved);
        else 
            newPatch = new Patch2D(params, smpi, ipatch, n_moved);
        
        // First, do a series of initializations
        newPatch->initStep1(params);
        newPatch->initStep2(params); // depends on 1D, 2D, ...
        newPatch->initStep3(params, smpi, n_moved);
        
        // clone vector of Species (virtual)
        newPatch->vecSpecies = SpeciesFactory::cloneVector(patch->vecSpecies, params, newPatch);
        
        // clone the electromagnetic fields (virtual)
        newPatch->EMfields   = ElectroMagnFactory::clone(patch->EMfields, params, newPatch->vecSpecies, newPatch);
        
        // interpolation operator (virtual)
        newPatch->Interp     = InterpolatorFactory::create(params, newPatch);
        // projection operator (virtual)
        newPatch->Proj       = ProjectorFactory::create(params, newPatch);
        
        // clone local diags
        newPatch->localDiags = DiagnosticFactory::cloneLocalDiagnostics(patch->localDiags, params, smpi, newPatch);
        newPatch->sio = SmileiIOFactory::clone(patch->sio, params, newPatch);
        
        // clone the collisions
        newPatch->vecCollisions = Collisions::clone(patch->vecCollisions, params);
        
        // clone the particle walls
        newPatch->partWalls = new PartWalls(patch->partWalls, newPatch);
        
        newPatch->createType(params);
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
