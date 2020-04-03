#ifndef PATCHESFACTORY_H
#define PATCHESFACTORY_H

#include "VectorPatch.h"
#include "Patch1D.h"
#include "Patch2D.h"
#include "Patch3D.h"
#include "PatchAM.h"
#include "DomainDecomposition.h"

#include "Tools.h"

#ifdef _VECTO
#include "SpeciesV.h"
#include "SpeciesVAdaptiveMixedSort.h"
#include "SpeciesVAdaptive.h"
#endif

class PatchesFactory
{
public:

    // Create one patch from scratch
    static Patch *create( Params &params, SmileiMPI *smpi, DomainDecomposition *domain_decomposition, unsigned int ipatch, unsigned int n_moved=0 )
    {
        if( params.geometry == "1Dcartesian" ) {
            return new Patch1D( params, smpi, domain_decomposition, ipatch, n_moved );
        } else if( params.geometry == "2Dcartesian" ) {
            return new Patch2D( params, smpi, domain_decomposition, ipatch, n_moved );
        } else if( params.geometry == "3Dcartesian" ) {
            return new Patch3D( params, smpi, domain_decomposition, ipatch, n_moved );
        } else if( params.geometry == "AMcylindrical" ) {
            return new PatchAM( params, smpi, domain_decomposition, ipatch, n_moved );
        }
        return nullptr;
    }
    
    // Clone one patch (avoid reading again the namelist)
    static Patch *clone( Patch *patch, Params &params, SmileiMPI *smpi, DomainDecomposition *domain_decomposition, unsigned int ipatch, unsigned int n_moved=0, bool with_particles = true )
    {
        if( params.geometry == "1Dcartesian" ) {
            return new Patch1D( static_cast<Patch1D *>( patch ), params, smpi, domain_decomposition, ipatch, n_moved, with_particles );
        } else if( params.geometry == "2Dcartesian" ) {
            return new Patch2D( static_cast<Patch2D *>( patch ), params, smpi, domain_decomposition, ipatch, n_moved, with_particles );
        } else if( params.geometry == "3Dcartesian" ) {
            return new Patch3D( static_cast<Patch3D *>( patch ), params, smpi, domain_decomposition, ipatch, n_moved, with_particles );
        } else if( params.geometry == "AMcylindrical" ) {
            return new PatchAM( static_cast<PatchAM *>( patch ), params, smpi, domain_decomposition, ipatch, n_moved, with_particles );
        }
        return nullptr;
    }
    
    // Create a vector of patches
    static void createVector( VectorPatch &vecPatches, Params &params, SmileiMPI *smpi, OpenPMDparams &openPMD, RadiationTables * radiation_tables_, unsigned int itime, unsigned int n_moved=0 )
    {
    
        vecPatches.diag_flag = ( params.restart? false : true );
        vecPatches.lastIterationPatchesMoved = itime;
        
        // Compute npatches (1 is std MPI behavior)
        unsigned int npatches, firstpatch;
        npatches = smpi->patch_count[smpi->getRank()];// Number of patches owned by current MPI process.
        firstpatch = 0;
        for( unsigned int impi = 0 ; impi < ( unsigned int )smpi->getRank() ; impi++ ) {
            firstpatch += smpi->patch_count[impi];
        }
        
        // If test mode, only 1 patch created
        if( smpi->test_mode ) {
            npatches = 1;
        }
        
        // Create patches (create patch#0 then clone it)
        vecPatches.resize( npatches );
        
        vecPatches.patches_[0] = create( params, smpi, vecPatches.domain_decomposition_, firstpatch, n_moved );
        
        TITLE( "Initializing Patches" );
        MESSAGE( 1, "First patch created" );
        
        // If normal mode (not test mode) clone the first patch to create the others
        unsigned int percent=10;
        for( unsigned int ipatch = 1 ; ipatch < npatches ; ipatch++ ) {
            if( ( 100*ipatch )/npatches > percent ) {
                MESSAGE( 2, "Approximately "<<percent<<"% of patches created" );
                percent += 10;
            }
            vecPatches.patches_[ipatch] = clone( vecPatches( 0 ), params, smpi, vecPatches.domain_decomposition_, firstpatch + ipatch, n_moved );
        }
        
        //Cleaning arrays and pointer
        for( unsigned int ispec=0 ; ispec<vecPatches( 0 )->vecSpecies.size(); ispec++ ) {
            //If a species was initialized via a numpy array
            if( vecPatches.patches_[0]->vecSpecies[ispec]->position_initialization_array_ ) {
                //delete the array
                delete vecPatches.patches_[0]->vecSpecies[ispec]->position_initialization_array_;
                //and never again create particles from this array. Pointer is kept to not NULL to remember this species was initialized from an array.
                for( unsigned int ipatch=0 ; ipatch < npatches ; ipatch++ ) {
                    vecPatches.patches_[ipatch]->vecSpecies[ispec]->n_numpy_particles_ = 0 ;
                }
            }
        }
        for( unsigned int ispec=0 ; ispec<vecPatches( 0 )->vecSpecies.size(); ispec++ ) {
            //If a species was initialized via a numpy array
            if( vecPatches.patches_[0]->vecSpecies[ispec]->momentum_initialization_array_ ) {
                //delete the array. Pointer is kept to not NULL to remember this species was initialized from an array.
                delete vecPatches.patches_[0]->vecSpecies[ispec]->momentum_initialization_array_;
            }
        }
        
        
        
        MESSAGE( 1, "All patches created" );
        
        vecPatches.setRefHindex();
        
        vecPatches.updateFieldList( smpi );
        
        TITLE( "Creating Diagnostics, antennas, and external fields" )
        vecPatches.createDiags( params, smpi, openPMD, radiation_tables_ );
        
        TITLE( "finalize MPI" )
        for( unsigned int ipatch = 0 ; ipatch < npatches ; ipatch++ ) {
            vecPatches.patches_[ipatch]->finalizeMPIenvironment( params );
        }
        vecPatches.nrequests = vecPatches( 0 )->requests_.size();
        
        
        // Figure out if there are antennas
        vecPatches.nAntennas = vecPatches( 0 )->EMfields->antennas.size();
        
        // Initialize lasers and antennas
        if( ! smpi->test_mode ) {
            vecPatches.initExternals( params );
        }
        
        MESSAGE( 1, "Done initializing diagnostics, antennas, and external fields" );
    }
    
};

#endif
