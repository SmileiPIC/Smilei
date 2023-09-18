#ifndef DIAGNOSTICFACTORY_H
#define DIAGNOSTICFACTORY_H

#include "VectorPatch.h"

#include "DiagnosticParticleBinning.h"
#include "DiagnosticScreen.h"
#include "DiagnosticRadiationSpectrum.h"
#include "DiagnosticProbes.h"
#include "DiagnosticScalar.h"
#include "DiagnosticTrack.h"
#include "DiagnosticNewParticles.h"
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

    static std::vector<Diagnostic *> createGlobalDiagnostics( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, RadiationTables * radiation_tables_ )
    {
        std::vector<Diagnostic *> vecDiagnostics;
        vecDiagnostics.push_back( new DiagnosticScalar( params, smpi, vecPatches( 0 ) ) );
        
        for( unsigned int i = 0, n = PyTools::nComponents( "DiagParticleBinning" ); i < n; i++ ) {
            vecDiagnostics.push_back( new DiagnosticParticleBinning( params, smpi, vecPatches( 0 ), i ) );
        }
        
        for( unsigned int i = 0, n = PyTools::nComponents( "DiagScreen" ); i < n; i++ ) {
            vecDiagnostics.push_back( new DiagnosticScreen( params, smpi, vecPatches( 0 ), i ) );
        }
        
        for( unsigned int i = 0, n = PyTools::nComponents( "DiagRadiationSpectrum" ); i < n; i++ ) {
            vecDiagnostics.push_back( new DiagnosticRadiationSpectrum(params, smpi, vecPatches(0), radiation_tables_ , i) );
        }
        
        return vecDiagnostics;
        
    } // END createGlobalDiagnostics
    
    
    
    static std::vector<Diagnostic *> createLocalDiagnostics( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, OpenPMDparams &openPMD )
    {
        std::vector<Diagnostic *> vecDiagnostics;
        
        for( unsigned int i = 0, n = PyTools::nComponents( "DiagFields" ); i < n; i++ ) {
            vecDiagnostics.push_back( DiagnosticFieldsFactory::create( params, smpi, vecPatches, i, openPMD ) );
        }
        
        for( unsigned int i = 0, n = PyTools::nComponents( "DiagProbe" ); i < n; i++ ) {
            vecDiagnostics.push_back( new DiagnosticProbes( params, smpi, vecPatches, i ) );
        }
        
        for( unsigned int i = 0, n = PyTools::nComponents( "DiagTrackParticles" ); i < n; i++ ) {
            vecDiagnostics.push_back( new DiagnosticTrack( params, smpi, vecPatches, i, vecDiagnostics.size(), openPMD ) );
        }
        
        for( unsigned int i = 0, n = PyTools::nComponents( "DiagNewParticles" ); i < n; i++ ) {
            vecDiagnostics.push_back( new DiagnosticNewParticles( params, smpi, vecPatches, i, vecDiagnostics.size(), openPMD ) );
        }
        
        if( PyTools::nComponents( "DiagPerformances" ) > 0 ) {
            vecDiagnostics.push_back( new DiagnosticPerformances( params, smpi ) );
        }
        
        return vecDiagnostics;
        
    } // END createLocalDiagnostics
    
    
    
    static std::vector<ProbeParticles *> createProbes()
    {
        std::vector<ProbeParticles *> probes( 0 );
        
        for( unsigned int i = 0, n = PyTools::nComponents( "DiagProbe" ); i < n; i++ ) {
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
