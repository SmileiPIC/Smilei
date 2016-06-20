#ifndef DIAGNOSTICFACTORY_H
#define DIAGNOSTICFACTORY_H

#include "DiagnosticParticles.h"
#include "DiagnosticProbes.h"
#include "DiagnosticScalar.h"
#include "DiagnosticTrack.h"
#include "DiagnosticTrackOld.h"

#include "DiagnosticFields1D.h"
#include "DiagnosticFields2D.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Create appropriate IO environment for the geometry 
//! \param params : Parameters
//! \param smpi : MPI environment
//  --------------------------------------------------------------------------------------------------------------------
class DiagnosticFieldsFactory {
public:
    static Diagnostic* create(Params& params, SmileiMPI* smpi, Patch* patch, unsigned int idiag) {
        Diagnostic* diag = NULL;
        if ( params.geometry == "1d3v" ) {
            diag = new DiagnosticFields1D(params, smpi, patch, idiag);
        }
        else if ( params.geometry == "2d3v" ) {
            diag = new DiagnosticFields2D(params, smpi, patch, idiag);
        }
        else {
            ERROR( "Geometry " << params.geometry << " not implemented" );
        }
        
        return diag;
    }
};

class DiagnosticFactory {
public:

    static std::vector<Diagnostic*> createGlobalDiagnostics(Params& params, SmileiMPI* smpi, Patch* patch) {
        std::vector<Diagnostic*> vecDiagnostics;
        vecDiagnostics.push_back( new DiagnosticScalar(params, smpi, patch, 0) ); // 0 : 1 scalar only (useless)
        
        for (unsigned int n_diag_particles = 0; n_diag_particles < PyTools::nComponents("DiagParticles"); n_diag_particles++) {
            vecDiagnostics.push_back( new DiagnosticParticles(params, smpi, patch, n_diag_particles) );
        }
        
        return vecDiagnostics;
    } // END createGlobalDiagnostics
    
    
    static std::vector<Diagnostic*> createLocalDiagnostics(Params& params, SmileiMPI* smpi, Patch* patch) {
        std::vector<Diagnostic*> vecDiagnostics;
        
        for (unsigned int n_diag_probes = 0; n_diag_probes < PyTools::nComponents("DiagProbe"); n_diag_probes++) {
            vecDiagnostics.push_back( new DiagnosticProbes(params, smpi, patch, vecDiagnostics.size(), n_diag_probes) );
        }
        
        return vecDiagnostics;
    } // END createLocalDiagnostics
    
    
    static std::vector<Diagnostic*> createOtherDiagnostics(Params& params, SmileiMPI* smpi, Patch* patch) {
        std::vector<Diagnostic*> vecDiagnostics;
        
        for (unsigned int n_diag_fields = 0; n_diag_fields < PyTools::nComponents("DiagFields"); n_diag_fields++) {
            vecDiagnostics.push_back( DiagnosticFieldsFactory::create(params, smpi, patch, n_diag_fields) );
        }
        
        for (unsigned int n_diag_track = 0; n_diag_track < patch->vecSpecies.size(); n_diag_track++) {
            if ( patch->vecSpecies[n_diag_track]->particles->tracked ) {
                vecDiagnostics.push_back( new DiagnosticTrack(params, smpi, patch, vecDiagnostics.size(), n_diag_track) );
            }
        }
        
        return vecDiagnostics;
    } // END createOtherDiagnostics
    
    
    // Cloning factory for local diags (global don't need cloning)
    static std::vector<Diagnostic*> cloneLocalDiagnostics(std::vector<Diagnostic*> vecDiagnostics, Params& params, SmileiMPI* smpi, Patch* patch) {
        std::vector<Diagnostic*> newVecDiagnostics(0);
        for (int idiag=0; idiag<vecDiagnostics.size(); idiag++) {
            if (vecDiagnostics[idiag]->type_ == "Probes" ) {
                newVecDiagnostics.push_back(
                    new DiagnosticProbes(static_cast<DiagnosticProbes*>(vecDiagnostics[idiag]), params, patch)
                );
            } else if (vecDiagnostics[idiag]->type_ == "Track" ) {
                newVecDiagnostics.push_back(
                    new DiagnosticTrackOld(static_cast<DiagnosticTrackOld*>(vecDiagnostics[idiag]), patch)
                );
            }
        }
        return newVecDiagnostics;
    }
    

};

#endif

