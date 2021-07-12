#ifndef VECTORPATCH_H
#define VECTORPATCH_H

#include <vector>
#include <iostream>
#include <cstdlib>
#include <iomanip>

#include "SpeciesFactory.h"
#include "InterpolatorFactory.h"
#include "ProjectorFactory.h"

#include "DiagnosticScalar.h"

#include "Checkpoint.h"
#include "OpenPMDparams.h"
#include "SmileiMPI.h"
#include "SimWindow.h"
#include "Timers.h"
#include "RadiationTables.h"
#include "ParticleCreator.h"

class Field;
class Timer;
class SimWindow;
class DomainDecomposition;

//! Class vectorPatch
//! This class corresponds to the MPI Patch Collection.

class VectorPatch
{
public :

    VectorPatch();
    VectorPatch( Params &params );
    ~VectorPatch();
    void saveOldRho( Params &params );
    void setMagneticFieldsForDiagnostic( Params &params );
    
    void close( SmileiMPI * );
    
    //! VectorPatch =
    //! - std::vector<Patch*>
    //! - interfaces between main programs & main PIC operators
    //! - methods to balance computation
    std::vector<Patch *> patches_;
    
    //! Vector of global diagnostics (diagnostics which cannot be computed locally)
    std::vector<Diagnostic *> globalDiags;
    //! Vector of local diagnostics (diagnostics which can partly be computed locally)
    std::vector<Diagnostic *> localDiags;
    
    //! Some vector operations extended to VectorPatch
    inline void resize( int npatches )
    {
        patches_.resize( npatches );
    }
    inline unsigned int  size() const
    {
        return patches_.size();
    }
    inline Patch *operator()( int ipatch )
    {
        return patches_[ipatch];
    }
    
    //! Set Id of the 1st patch stored on the current MPI process
    //!   used during balancing
    inline void setRefHindex()
    {
        refHindex_ = patches_[0]->Hindex();
    }
    //! Resize vector of field*
    void updateFieldList( SmileiMPI *smpi );
    void updateFieldList( int ispec, SmileiMPI *smpi );
    
    //! Create the diagnostic list
    void createDiags( Params &params, SmileiMPI *smpi, OpenPMDparams &, RadiationTables * radiation_tables_ );
    
    //! get a particular scalar
    inline double getScalar( std::string name )
    {
        DiagnosticScalar *diag = static_cast<DiagnosticScalar *>( globalDiags[0] );
        return diag->getScalar( name );
    }
    
    bool needsRhoJsNow( int timestep )
    {
        // Figure out whether scalars need Rho and Js
        if( globalDiags[0]->needsRhoJs( timestep ) ) {
            return true;
        }
        
        // Figure out whether fields or probes need Rho and Js
        for( unsigned int i=0; i<localDiags.size(); i++ )
            if( localDiags[i]->needsRhoJs( timestep ) ) {
                return true;
            }
            
        return false;
    }
    
    // Interfaces between main programs & main PIC operators
    // -----------------------------------------------------
    
    //! Reconfigure all patches for the new time step
    void configuration( Params &params, Timers &timers, int itime );
    
    //! Reconfigure all patches for the new time step
    void reconfiguration( Params &params, Timers &timers, int itime );
    
    //! Particle sorting for all patches
    void sortAllParticles( Params &params );
    
    //! For all patch, move particles (restartRhoJ(s), dynamics and exchangeParticles)
    void dynamics( Params &params,
                   SmileiMPI *smpi,
                   SimWindow *simWindow,
                   RadiationTables &RadiationTables,
                   MultiphotonBreitWheelerTables &MultiphotonBreitWheelerTables,
                   double time_dual,
                   Timers &timers, int itime );
    
    //! For all patches, exchange particles and sort them.
    void finalizeAndSortParticles( Params &params, SmileiMPI *smpi, SimWindow *simWindow,
                                  double time_dual,
                                  Timers &timers, int itime );
    void finalizeSyncAndBCFields( Params &params, SmileiMPI *smpi, SimWindow *simWindow,
                                      double time_dual, Timers &timers, int itime );

    //! Particle merging
    void mergeParticles(Params &params, SmileiMPI *smpi, double time_dual,Timers &timers, int itime );

    //! Clean MPI buffers and resize particle arrays to save memory
    void cleanParticlesOverhead(Params &params, Timers &timers, int itime );
                              
    //! Particle injection from the boundaries
    void injectParticlesFromBoundaries( Params &params, Timers &timers, unsigned int itime );
                                      
    //! Computation of the total charge
    void computeCharge(bool old = false);
    
    void projectionForDiags( Params &params,
                               SmileiMPI *smpi,
                               SimWindow *simWindow,
                               double time_dual,
                               Timers &timers, int itime );
                               
    // compute rho only given by relativistic species which require initialization of the relativistic fields
    void computeChargeRelativisticSpecies( double time_primal, Params &params );
    
    // run particles ponderomptive dynamics, envelope's solver
    void runEnvelopeModule( Params &params,
            SmileiMPI *smpi,
            SimWindow *simWindow,
            double time_dual, Timers &timers, int itime );
    //! For all patches, deposit susceptibility, then advance momentum of particles interacting with envelope
    void ponderomotiveUpdateSusceptibilityAndMomentum( Params &params,
            SmileiMPI *smpi,
            SimWindow *simWindow,
            double time_dual, Timers &timers, int itime );
    //! For all patches, advance position of particles interacting with envelope, comm particles, project charge and current density
    void ponderomotiveUpdatePositionAndCurrents( Params &params,
            SmileiMPI *smpi,
            SimWindow *simWindow,
            double time_dual, Timers &timers, int itime );
    void resetRhoJ(bool old = false);
    
    //! For all patch, sum densities on ghost cells (sum per species if needed, sync per patch and MPI sync)
    void sumDensities( Params &params, double time_dual, Timers &timers, int itime, SimWindow *simWindow, SmileiMPI *smpi );
    
    //! For all patch, sum susceptibility on ghost cells (sum per species if needed, sync per patch and MPI sync)
    void sumSusceptibility( Params &params, double time_dual, Timers &timers, int itime, SimWindow *simWindow, SmileiMPI *smpi );
    
    //! For all patch, update E and B (Ampere, Faraday, boundary conditions, exchange B and center B)
    void solveMaxwell( Params &params, SimWindow *simWindow, int itime, double time_dual,
                       Timers &timers, SmileiMPI *smpi );
                       
    //! For all patch, update envelope field A (envelope equation, boundary contitions, exchange A)
    void solveEnvelope( Params &params, SimWindow *simWindow, int itime, double time_dual, Timers &timers, SmileiMPI *smpi );
    
    //! For all patch, Compute and Write all diags (Scalars, Probes, Phases, TrackParticles, Fields, Average fields)
    void runAllDiags( Params &params, SmileiMPI *smpi, unsigned int itime, Timers &timers, SimWindow *simWindow );
    void initAllDiags( Params &params, SmileiMPI *smpi );
    void closeAllDiags( SmileiMPI *smpi );
    
    //! Check if rho is null (MPI & patch sync)
    bool isRhoNull( SmileiMPI *smpi );
    
    //! Solve Poisson to initialize E
    void solvePoisson( Params &params, SmileiMPI *smpi );
    void runNonRelativisticPoissonModule( Params &params, SmileiMPI* smpi,  Timers &timers );
    void solvePoissonAM( Params &params, SmileiMPI *smpi);
    
    //! Solve relativistic Poisson problem to initialize E and B of a relativistic bunch
    void runRelativisticModule( double time_prim, Params &params, SmileiMPI* smpi,  Timers &timers );
    void solveRelativisticPoisson( Params &params, SmileiMPI *smpi, double time_primal );
    void solveRelativisticPoissonAM( Params &params, SmileiMPI *smpi, double time_primal );
    
    //! For all patch initialize the externals (lasers, fields, antennas)
    void initExternals( Params &params );
    
    //! For all patches, apply the antenna current
    void applyAntennas( double time );
    
    //! For all patches, apply collisions
    void applyCollisions( Params &params, int itime, Timers &timer );
    
    //! For all patches, allocate a field if not allocated
    void allocateField( unsigned int ifield, Params &params );
    
    //! For each patch, apply external fields
    void applyExternalFields();
    
    //! For each patch, apply external time fields
    void applyPrescribedFields(double time);

	//! reset all external time fields;
    void resetPrescribedFields();
    
    void saveExternalFields( Params &params );
    
    //  Balancing methods
    // ------------------
    
    //! Wrapper of load balancing methods, including SmileiMPI::recompute_patch_count. Called from main program
    void loadBalance( Params &params, double time_dual, SmileiMPI *smpi, SimWindow *simWindow, unsigned int itime );
    
    //! Explicits patch movement regarding new patch distribution stored in smpi->patch_count
    void createPatches( Params &params, SmileiMPI *smpi, SimWindow *simWindow );
    
    //! Exchange patches, based on createPatches initialization
    void exchangePatches( SmileiMPI *smpi, Params &params );
    
    //! Write in a file patches communications
    void outputExchanges( SmileiMPI *smpi );
    
    //! Init new envelope from input namelist
    void initNewEnvelope( Params &params );
    
    // Lists of fields
    std::vector<Field *> densities;
    
    std::vector<Field *> Bs0;
    std::vector<Field *> Bs1;
    std::vector<Field *> Bs2;
    std::vector<Field *> densitiesLocalx;
    std::vector<Field *> densitiesLocaly;
    std::vector<Field *> densitiesLocalz;
    std::vector<Field *> densitiesMPIx;
    std::vector<Field *> densitiesMPIy;
    std::vector<Field *> densitiesMPIz;
    
    std::vector<int> LocalxIdx;
    std::vector<int> LocalyIdx;
    std::vector<int> LocalzIdx;
    std::vector<int> MPIxIdx;
    std::vector<int> MPIyIdx;
    std::vector<int> MPIzIdx;
    
    std::vector<Field *> B_localx;
    std::vector<Field *> B_MPIx;
    
    std::vector<Field *> B1_localy;
    std::vector<Field *> B1_MPIy;
    
    std::vector<Field *> B2_localz;
    std::vector<Field *> B2_MPIz;
    
    std::vector<Field *> listJx_;
    std::vector<Field *> listJy_;
    std::vector<Field *> listJz_;
    std::vector<Field *> listrho_;
    std::vector<Field *> listJxs_;
    std::vector<Field *> listJys_;
    std::vector<Field *> listJzs_;
    std::vector<Field *> listrhos_;
    std::vector<Field *> listEx_;
    std::vector<Field *> listEy_;
    std::vector<Field *> listEz_;
    std::vector<Field *> listBx_;
    std::vector<Field *> listBy_;
    std::vector<Field *> listBz_;
    
    std::vector<Field *> listA_;
    std::vector<Field *> listA0_;
    // std::vector<Field *> listEnvE_;
    std::vector<Field *> listEnvEx_;
    // std::vector<Field *> listEnvA_;
    // std::vector<Field *> listPhi_;
    // std::vector<Field *> listPhi0_;
    std::vector<Field *> listGradPhix_;
    std::vector<Field *> listGradPhiy_;
    std::vector<Field *> listGradPhiz_;
    std::vector<Field *> listGradPhil_;
    std::vector<Field *> listGradPhir_;
    std::vector<Field *> listGradPhix0_;
    std::vector<Field *> listGradPhiy0_;
    std::vector<Field *> listGradPhiz0_;
    std::vector<Field *> listGradPhil0_;
    std::vector<Field *> listGradPhir0_;
    std::vector<Field *> listEnv_Chi_;
    std::vector<Field *> listEnv_Chis_;
    
    std::vector<std::vector< Field *>> listJl_;
    std::vector<std::vector< Field *>> listJr_;
    std::vector<std::vector< Field *>> listJt_;
    std::vector<std::vector< Field *>> listrho_AM_;
    std::vector<std::vector< Field *>> listrho_old_AM_;
    std::vector<std::vector< Field *>> listJls_;
    std::vector<std::vector< Field *>> listJrs_;
    std::vector<std::vector< Field *>> listJts_;
    std::vector<std::vector< Field *>> listrhos_AM_;
    std::vector<std::vector< Field *>> listEl_;
    std::vector<std::vector< Field *>> listEr_;
    std::vector<std::vector< Field *>> listEt_;
    std::vector<std::vector< Field *>> listBl_;
    std::vector<std::vector< Field *>> listBr_;
    std::vector<std::vector< Field *>> listBt_;
    
    
    //! True if any antennas
    unsigned int nAntennas;
    
    //! 1st patch index of patches_ (stored for balancing op)
    int refHindex_;
    
    //! Count global (MPI x patches) number of particles
    uint64_t getGlobalNumberOfParticles( SmileiMPI *smpi )
    {
        std::vector<uint64_t> nParticles = getGlobalNumberOfParticlesPerSpecies( smpi );
        
        uint64_t global_number_particles = 0;
        for( unsigned int i = 0; i<nParticles.size(); i++ ) {
            global_number_particles += nParticles[i];
        }
        return global_number_particles;
    }
    //! Count global (MPI x patches) number of particles per species
    std::vector<uint64_t> getGlobalNumberOfParticlesPerSpecies( SmileiMPI *smpi )
    {
        unsigned int nSpecies( ( *this )( 0 )->vecSpecies.size() );
        std::vector< uint64_t > nParticles( nSpecies, 0 );
        for( unsigned int ipatch = 0 ; ipatch < this->size() ; ipatch++ ) {
            for( unsigned int ispec = 0 ; ispec < nSpecies ; ispec++ ) {
                nParticles[ispec] += ( *this )( ipatch )->vecSpecies[ispec]->getNbrOfParticles();
            }
        }
        for( unsigned int ispec = 0 ; ispec < nSpecies ; ispec++ ) {
            uint64_t tmp( 0 );
            MPI_Reduce( &( nParticles[ispec] ), &tmp, 1, MPI_UINT64_T, MPI_SUM, 0, smpi->world() );
            nParticles[ispec] = tmp;
        }
        return nParticles;
    }
    //! Print global (MPI x patches) number of particles per species
    void printGlobalNumberOfParticlesPerSpecies( SmileiMPI *smpi )
    {
        std::vector< uint64_t > nParticles = getGlobalNumberOfParticlesPerSpecies( smpi );
        for( unsigned int ispec = 0 ; ispec < nParticles.size() ; ispec++ ) {
            MESSAGE( 2, "Species " << ispec << " (" << ( *this )( 0 )->vecSpecies[ispec]->name_ << ") created with " << nParticles[ispec] << " particles" );
        }
    }
    
    void checkMemoryConsumption( SmileiMPI *smpi, VectorPatch *region_vecpatches );
    
    void checkExpectedDiskUsage( SmileiMPI *smpi, Params &params, Checkpoint &checkpoint );
    
    // Keep track if we need the needsRhoJsNow
    int diag_flag;
    
    int nrequests;
    
    //! Tells which iteration was last time the patches moved (by moving window or load balancing)
    unsigned int lastIterationPatchesMoved;
    
    DomainDecomposition *domain_decomposition_;
    
    
    //! Methods to access readably to patch PIC operators.
    //!   - patches_ should not be access outsied of VectorPatch
    //!   - for now in SimWindow
    inline Species *species( int ipatch, int ispec )
    {
        return ( *this )( ipatch )->vecSpecies[ispec];
    }
    
    inline Interpolator *inter( int ipatch, int ispec )
    {
        return ( *this )( ipatch )->vecSpecies[ispec]->Interp;
    }
    
    inline ElectroMagn *emfields( int ipatch )
    {
        return ( *this )( ipatch )->EMfields;
    }
    
    inline Projector *proj( int ipatch, int ispec )
    {
        return ( *this )( ipatch )->vecSpecies[ispec]->Proj;
    }
    
    inline PartWalls *partwalls( int ipatch )
    {
        return ( *this )( ipatch )->partWalls;
    }
    
private :

    //  Internal balancing members
    // ---------------------------
    std::vector<Patch *> recv_patches_;
    
    std::vector<int> recv_patch_id_;
    std::vector<int> send_patch_id_;
    
    //! Current intensity of antennas
    double antenna_intensity;
    
    std::vector<Timer *> diag_timers;
};


#endif
