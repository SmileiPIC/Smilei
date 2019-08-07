#ifndef DIAGNOSTICFACTORY_H
#define DIAGNOSTICFACTORY_H

#include "VectorPatch.h"

#include "DiagnosticParticleBinning.h"
#include "DiagnosticScreen.h"
#include "DiagnosticRadiationSpectrum.h"
#include "DiagnosticProbes.h"
#include "DiagnosticScalar.h"
#include "DiagnosticTrack.h"
#include "DiagnosticPerformances.h"

#include "DiagnosticFields1D.h"
#include "DiagnosticFields2D.h"
#include "DiagnosticFields3D.h"
#include "DiagnosticFieldsAM.h"

//  --------------------------------------------------------------------------------------------------------------------
//! Create appropriate IO environment for the geometry
//! \param params : Parameters
//! \param smpi : MPI environment
//  --------------------------------------------------------------------------------------------------------------------
class DiagnosticFieldsFactory
{
public:
    static Diagnostic *create( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, unsigned int idiag, OpenPMDparams &openPMD )
    {
        Diagnostic *diag = NULL;
        if( params.geometry == "1Dcartesian" ) {
            diag = new DiagnosticFields1D( params, smpi, vecPatches, idiag, openPMD );
        } else if( params.geometry == "2Dcartesian" ) {
            diag = new DiagnosticFields2D( params, smpi, vecPatches, idiag, openPMD );
        } else if( params.geometry == "3Dcartesian" ) {
            diag = new DiagnosticFields3D( params, smpi, vecPatches, idiag, openPMD );
        } else if( params.geometry == "AMcylindrical" ) {
            diag = new DiagnosticFieldsAM( params, smpi, vecPatches, idiag, openPMD );
        } else {
            ERROR( "Geometry " << params.geometry << " not implemented" );
        }
        
        return diag;
    }
};

class DiagnosticFactory
{
public:

    static std::vector<Diagnostic *> createGlobalDiagnostics( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches )
    {
        std::vector<Diagnostic *> vecDiagnostics;
        vecDiagnostics.push_back( new DiagnosticScalar( params, smpi, vecPatches( 0 ) ) );
        
        for( unsigned int n_diag_particles = 0; n_diag_particles < PyTools::nComponents( "DiagParticleBinning" ); n_diag_particles++ ) {
            vecDiagnostics.push_back( new DiagnosticParticleBinning( params, smpi, vecPatches( 0 ), n_diag_particles ) );
        }
        
        for( unsigned int n_diag_screen = 0; n_diag_screen < PyTools::nComponents( "DiagScreen" ); n_diag_screen++ ) {
            vecDiagnostics.push_back( new DiagnosticScreen( params, smpi, vecPatches( 0 ), n_diag_screen ) );
        }

        for (unsigned int n_diag_rad_spectrum = 0; n_diag_rad_spectrum < PyTools::nComponents("DiagRadiationSpectrum"); n_diag_rad_spectrum++) {
            vecDiagnostics.push_back( new DiagnosticRadiationSpectrum(params, smpi, vecPatches(0), n_diag_rad_spectrum) );
        }
        
//MESSAGE ("Glob diag");
        return vecDiagnostics;
        
    } // END createGlobalDiagnostics
    
    
    
    static std::vector<Diagnostic *> createLocalDiagnostics( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, OpenPMDparams &openPMD )
    {
        std::vector<Diagnostic *> vecDiagnostics;
        //MESSAGE("in create local diags:  global dims after declaring vecdiag " << vecPatches(0)->EMfields->Jx_s[1]->globalDims_);
        
        for( unsigned int n_diag_fields = 0; n_diag_fields < PyTools::nComponents( "DiagFields" ); n_diag_fields++ ) {
            vecDiagnostics.push_back( DiagnosticFieldsFactory::create( params, smpi, vecPatches, n_diag_fields, openPMD ) );
            //  MESSAGE("in create local diags:  global dims after creating and pushing back field diag " << vecPatches(0)->EMfields->Jx_s[1]->globalDims_);
        }
        
        for( unsigned int n_diag_probe = 0; n_diag_probe < PyTools::nComponents( "DiagProbe" ); n_diag_probe++ ) {
            vecDiagnostics.push_back( new DiagnosticProbes( params, smpi, vecPatches, n_diag_probe ) );
        }
        
        for( unsigned int n_diag_track = 0; n_diag_track < PyTools::nComponents( "DiagTrackParticles" ); n_diag_track++ ) {
            vecDiagnostics.push_back( new DiagnosticTrack( params, smpi, vecPatches, n_diag_track, vecDiagnostics.size(), openPMD ) );
        }
        
        if( PyTools::nComponents( "DiagPerformances" ) > 0 ) {
            vecDiagnostics.push_back( new DiagnosticPerformances( params, smpi ) );
        }
        
//MESSAGE("local diag");
        return vecDiagnostics;
        
    } // END createLocalDiagnostics
    
    
    
    static std::vector<ProbeParticles *> createProbes()
    {
        std::vector<ProbeParticles *> probes( 0 );
        
        for( unsigned int n_probe = 0; n_probe < PyTools::nComponents( "DiagProbe" ); n_probe++ ) {
            probes.push_back( new ProbeParticles() );
        }
        
        return probes;
    } // END createProbes
    
    
    static std::vector<ProbeParticles *> cloneProbes( std::vector<ProbeParticles *> probes )
    {
        std::vector<ProbeParticles *> newProbes( 0 );
        
        for( unsigned int n_probe=0; n_probe<probes.size(); n_probe++ ) {
            newProbes.push_back( new ProbeParticles( probes[n_probe] ) );
        }
        
        return newProbes;
    }
    
    
};

#endif
