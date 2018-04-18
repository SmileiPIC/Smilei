#ifndef PATCHESFACTORY_H
#define PATCHESFACTORY_H

#include "VectorPatch.h"
#include "Patch1D.h"
#include "Patch2D.h"
#include "Patch3D.h"
#include "DomainDecomposition.h"

#include "Tools.h"

#ifdef _VECTO
#include "SpeciesV.h"
#include "SpeciesDynamicV.h"
#include "SpeciesDynamicV2.h"
#endif

class PatchesFactory {
public:

    // Create one patch from scratch
    static Patch* create(Params& params, SmileiMPI* smpi, DomainDecomposition* domain_decomposition, unsigned int ipatch, unsigned int n_moved=0) {
        if (params.geometry == "1Dcartesian")
            return new Patch1D(params, smpi, domain_decomposition, ipatch, n_moved);
        else if (params.geometry == "2Dcartesian")
            return new Patch2D(params, smpi, domain_decomposition, ipatch, n_moved);
        else if (params.geometry == "3Dcartesian")
            return new Patch3D(params, smpi, domain_decomposition, ipatch, n_moved);
        return nullptr;
    }

    // Clone one patch (avoid reading again the namelist)
    static Patch* clone(Patch* patch, Params& params, SmileiMPI* smpi, DomainDecomposition* domain_decomposition, unsigned int ipatch, unsigned int n_moved=0, bool with_particles = true) {
        if (params.geometry == "1Dcartesian")
            return new Patch1D(static_cast<Patch1D*>(patch), params, smpi, domain_decomposition, ipatch, n_moved, with_particles);
        else if (params.geometry == "2Dcartesian")
            return new Patch2D(static_cast<Patch2D*>(patch), params, smpi, domain_decomposition, ipatch, n_moved, with_particles);
        else if (params.geometry == "3Dcartesian")
            return new Patch3D(static_cast<Patch3D*>(patch), params, smpi, domain_decomposition, ipatch, n_moved, with_particles);
        return nullptr;
    }

    // Create a vector of patches
    static VectorPatch createVector(Params& params, SmileiMPI* smpi, OpenPMDparams& openPMD, unsigned int itime, unsigned int n_moved=0) {
        VectorPatch vecPatches( params );

        vecPatches.diag_flag = (params.restart? false : true);
        vecPatches.lastIterationPatchesMoved = itime;

        // Compute npatches (1 is std MPI behavior)
        unsigned int npatches, firstpatch;
        npatches = smpi->patch_count[smpi->getRank()];// Number of patches owned by current MPI process.
        firstpatch = 0;
        for (unsigned int impi = 0 ; impi < (unsigned int)smpi->getRank() ; impi++) {
            firstpatch += smpi->patch_count[impi];
        }
        DEBUG( smpi->getRank() << ", nPatch = " << npatches << " - starting at " << firstpatch );

        // If test mode, only 1 patch created
        if( smpi->test_mode ) npatches = 1;

        // Create patches (create patch#0 then clone it)
        vecPatches.resize(npatches);
        vecPatches.patches_[0] = create(params, smpi, vecPatches.domain_decomposition_, firstpatch, n_moved);

        TITLE("Initializing Patches");
        MESSAGE(1,"First patch created");

        // If normal mode (not test mode) clone the first patch to create the others
        unsigned int percent=10;
        for (unsigned int ipatch = 1 ; ipatch < npatches ; ipatch++) {
            if( (100*ipatch)/npatches > percent ) {
                MESSAGE(2,"Approximately "<<percent<<"% of patches created");
                percent += 10;
            }
            vecPatches.patches_[ipatch] = clone(vecPatches(0), params, smpi, vecPatches.domain_decomposition_, firstpatch + ipatch, n_moved);
        }

#ifdef _VECTO
        if (params.vecto == "normal") {
            //Need to sort because particles are not well sorted at creation
            for (unsigned int ipatch=0 ; ipatch < npatches ; ipatch++){
                for (unsigned int ispec=0 ; ispec<vecPatches(ipatch)->vecSpecies.size(); ispec++) {
                    if ( dynamic_cast<SpeciesV*>(vecPatches.patches_[ipatch]->vecSpecies[ispec]) )
                        dynamic_cast<SpeciesV*>(vecPatches.patches_[ipatch]->vecSpecies[ispec])->compute_part_cell_keys(params);
                    vecPatches.patches_[ipatch]->vecSpecies[ispec]->sort_part(params);
                }
            }
        }
        else if (params.vecto == "dynamic") {
            //Need to sort because particles are not well sorted at creation
            for (unsigned int ipatch=0 ; ipatch < npatches ; ipatch++){
                for (unsigned int ispec=0 ; ispec<vecPatches(ipatch)->vecSpecies.size(); ispec++) {
                    if ( dynamic_cast<SpeciesDynamicV*>(vecPatches.patches_[ipatch]->vecSpecies[ispec]) )
                    {
                        dynamic_cast<SpeciesDynamicV*>(vecPatches.patches_[ipatch]->vecSpecies[ispec])->compute_part_cell_keys(params);
                    }
                    if (dynamic_cast<SpeciesDynamicV*>(vecPatches.patches_[ipatch]->vecSpecies[ispec])->vectorized_operators)
                    {
                        vecPatches.patches_[ipatch]->vecSpecies[ispec]->sort_part(params);
                    }
                    else
                    {
                        vecPatches.patches_[ipatch]->vecSpecies[ispec]->Species::sort_part(params);
                    }
                }
            }
        }
        else if (params.vecto == "dynamic2") {
            //Need to sort because particles are not well sorted at creation
            for (unsigned int ipatch=0 ; ipatch < npatches ; ipatch++){
                for (unsigned int ispec=0 ; ispec<vecPatches(ipatch)->vecSpecies.size(); ispec++) {
                    if ( dynamic_cast<SpeciesDynamicV2*>(vecPatches.patches_[ipatch]->vecSpecies[ispec]) )
                    {
                        dynamic_cast<SpeciesDynamicV2*>(vecPatches.patches_[ipatch]->vecSpecies[ispec])->compute_part_cell_keys(params);
                    }
                    dynamic_cast<SpeciesDynamicV2*>(vecPatches.patches_[ipatch]->vecSpecies[ispec])->sort_part(params);
                }
            }
        }
#endif


        MESSAGE(1,"All patches created");

        vecPatches.set_refHindex();

        vecPatches.update_field_list();

        TITLE("Creating Diagnostics, antennas, and external fields")
        vecPatches.createDiags( params, smpi, openPMD );

        for (unsigned int ipatch = 0 ; ipatch < npatches ; ipatch++)
            vecPatches.patches_[ipatch]->finalizeMPIenvironment(params);
        vecPatches.nrequests = vecPatches(0)->requests_.size();


        // Figure out if there are antennas
        vecPatches.nAntennas = vecPatches(0)->EMfields->antennas.size();
        vecPatches.initExternals( params );

        MESSAGE(1,"Done initializing diagnostics, antennas, and external fields");
        return vecPatches;
    }

};

#endif
