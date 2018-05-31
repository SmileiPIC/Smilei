#ifndef VECTORPATCH_H
#define VECTORPATCH_H

#include <vector>
#include <iostream>
#include <cstdlib>
#include <iomanip>

#include "SpeciesFactory.h"
#include "ElectroMagnFactory.h"
#include "InterpolatorFactory.h"
#include "ProjectorFactory.h"

#include "DiagnosticScalar.h"

#include "Checkpoint.h"
#include "OpenPMDparams.h"
#include "SmileiMPI.h"
#include "SimWindow.h"
#include "Timers.h"
#include "RadiationTables.h"

class Field;
class Timer;
class SimWindow; 
class DomainDecomposition;

//! Class Patch : sub MPI domain
//!     Collection of patch = MPI domain

class VectorPatch {
public :

    VectorPatch();
    VectorPatch( Params &params );
    ~VectorPatch();
    void save_old_rho(Params &params); 

    void close(SmileiMPI*);

    //! VectorPatch =
    //! - std::vector<Patch*>
    //! - interfaces between main programs & main PIC operators
    //! - methods to balance computation
    std::vector<Patch*> patches_;

    //! Vector of global diagnostics (diagnostics which cannot be computed locally)
    std::vector<Diagnostic*> globalDiags;
    //! Vector of local diagnostics (diagnostics which can partly be computed locally)
    std::vector<Diagnostic*> localDiags;

    //! Some vector operations extended to VectorPatch
    inline void resize(int npatches) {
        patches_.resize(npatches);
    }
    inline unsigned int  size() const {
        return patches_.size();
    }
    inline Patch* operator()(int ipatch) {
        return patches_[ipatch];
    }

    //! Set Id of the 1st patch stored on the current MPI process
    //!   used during balancing
    inline void set_refHindex() {
        refHindex_ = patches_[0]->Hindex();
    }
    //! Resize vector of field*
    void update_field_list();
    void update_field_list(int ispec);

    void createDiags(Params& params, SmileiMPI* smpi, OpenPMDparams&);

    //! get a particular scalar
    inline double getScalar(std::string name) {
        DiagnosticScalar* diag = static_cast<DiagnosticScalar*>( globalDiags[0] );
        return diag->getScalar( name );
    }

    bool needsRhoJsNow( int timestep ) {
        // Figure out whether scalars need Rho and Js
        if( globalDiags[0]->needsRhoJs(timestep) ) return true;

        // Figure out whether fields or probes need Rho and Js
        for( unsigned int i=0; i<localDiags.size(); i++ )
            if( localDiags[i]->needsRhoJs(timestep) ) return true;

        return false;
    }

    // Interfaces between main programs & main PIC operators
    // -----------------------------------------------------

    //! For all patch, move particles (restartRhoJ(s), dynamics and exchangeParticles)
    void dynamics(Params& params,
                  SmileiMPI* smpi,
                  SimWindow* simWindow,
                  RadiationTables & RadiationTables,
                  MultiphotonBreitWheelerTables & MultiphotonBreitWheelerTables,
                  double time_dual,
                  Timers &timers, int itime);

    void finalize_and_sort_parts(Params& params, SmileiMPI* smpi, SimWindow* simWindow, 
                  RadiationTables & RadiationTables,
                  MultiphotonBreitWheelerTables & MultiphotonBreitWheelerTables,
                  double time_dual,
                  Timers &timers, int itime);
    void finalize_sync_and_bc_fields(Params& params, SmileiMPI* smpi, SimWindow* simWindow,
                  double time_dual, Timers &timers, int itime);

    void computeCharge();

    void projection_for_diags(Params& params,
                  SmileiMPI* smpi,
                  SimWindow* simWindow,
                  double time_dual,
                  Timers &timers, int itime);
 
    // compute rho only given by relativistic species which require initialization of the relativistic fields
    void computeChargeRelativisticSpecies(double time_primal); 

    void resetRhoJ();

    //! For all patch, sum densities on ghost cells (sum per species if needed, sync per patch and MPI sync)
    void sumDensities(Params &params, double time_dual, Timers &timers, int itime, SimWindow* simWindow );

    //! For all patch, update E and B (Ampere, Faraday, boundary conditions, exchange B and center B)
    void solveMaxwell(Params& params, SimWindow* simWindow, int itime, double time_dual,
                      Timers & timers);
    
    //! For all patch, Compute and Write all diags (Scalars, Probes, Phases, TrackParticles, Fields, Average fields)
    void runAllDiags(Params& params, SmileiMPI* smpi, unsigned int itime, Timers & timers, SimWindow* simWindow);
    void initAllDiags(Params& params, SmileiMPI* smpi);
    void closeAllDiags(SmileiMPI* smpi);
    void openAllDiags(Params& params, SmileiMPI* smpi);

    //! Check if rho is null (MPI & patch sync)
    bool isRhoNull( SmileiMPI* smpi );

    //! Solve Poisson to initialize E
    void solvePoisson( Params &params, SmileiMPI* smpi );

    //! Solve relativistic Poisson problem to initialize E and B of a relativistic bunch
    void solveRelativisticPoisson( Params &params, SmileiMPI* smpi, double time_primal );

    //! For all patch initialize the externals (lasers, fields, antennas)
    void initExternals(Params& params);

    //! For all patches, apply the antenna current
    void applyAntennas(double time);

    //! For all patches, apply collisions
    void applyCollisions(Params &params, int itime, Timers & timer);

    //! For each patch, apply external fields
    void applyExternalFields();

    void saveExternalFields( Params &params );

    //  Balancing methods
    // ------------------

    //! Wrapper of load balancing methods, including SmileiMPI::recompute_patch_count. Called from main program
    void load_balance(Params& params, double time_dual, SmileiMPI* smpi, SimWindow* simWindow, unsigned int itime);

    //! Explicits patch movement regarding new patch distribution stored in smpi->patch_count
    void createPatches(Params& params, SmileiMPI* smpi, SimWindow* simWindow);

    //! Exchange patches, based on createPatches initialization
    void exchangePatches(SmileiMPI* smpi, Params& params);

    //! Write in a file patches communications
    void output_exchanges(SmileiMPI* smpi);

    // Lists of fields
    std::vector<Field*> densities;

    std::vector<Field*> Bs0;
    std::vector<Field*> Bs1;
    std::vector<Field*> Bs2;
    std::vector<Field*> densitiesLocalx;
    std::vector<Field*> densitiesLocaly;
    std::vector<Field*> densitiesLocalz;
    std::vector<Field*> densitiesMPIx;
    std::vector<Field*> densitiesMPIy;
    std::vector<Field*> densitiesMPIz;

    std::vector<int> LocalxIdx;
    std::vector<int> LocalyIdx;
    std::vector<int> LocalzIdx;
    std::vector<int> MPIxIdx;
    std::vector<int> MPIyIdx;
    std::vector<int> MPIzIdx;

    std::vector<Field*> B_localx;
    std::vector<Field*> B_MPIx;

    std::vector<Field*> B1_localy;
    std::vector<Field*> B1_MPIy;

    std::vector<Field*> B2_localz;
    std::vector<Field*> B2_MPIz;

    std::vector<Field*> listJx_;
    std::vector<Field*> listJy_;
    std::vector<Field*> listJz_;
    std::vector<Field*> listrho_;
    std::vector<Field*> listJxs_;
    std::vector<Field*> listJys_;
    std::vector<Field*> listJzs_;
    std::vector<Field*> listrhos_;
    std::vector<Field*> listEx_;
    std::vector<Field*> listEy_;
    std::vector<Field*> listEz_;
    std::vector<Field*> listBx_;
    std::vector<Field*> listBy_;
    std::vector<Field*> listBz_;

    //! True if any antennas
    unsigned int nAntennas;

    //! 1st patch index of patches_ (stored for balancing op)
    int refHindex_;

    //! Count global (MPI x patches) number of particles per species
    void printNumberOfParticles(SmileiMPI* smpi) {
        unsigned int nSpecies( (*this)(0)->vecSpecies.size() );
        std::vector< uint64_t > nParticles( nSpecies, 0 );
        for (unsigned int ipatch = 0 ; ipatch < this->size() ; ipatch++ ) {
            for (unsigned int ispec = 0 ; ispec < nSpecies ; ispec++ ) {
                nParticles[ispec] += (*this)(ipatch)->vecSpecies[ispec]->getNbrOfParticles();
            }
        }
        for (unsigned int ispec = 0 ; ispec < nSpecies ; ispec++ ) {
            uint64_t tmp(0);
            MPI_Reduce( &(nParticles[ispec]), &tmp, 1, MPI_UINT64_T , MPI_SUM, 0, smpi->SMILEI_COMM_WORLD );
            MESSAGE(2, "Species " << ispec << " (" << (*this)(0)->vecSpecies[ispec]->name << ") created with " << tmp << " particles" );
        }
    }

    void check_memory_consumption(SmileiMPI* smpi);
    
    void check_expected_disk_usage( SmileiMPI* smpi, Params& params, Checkpoint& checkpoint);

    // Keep track if we need the needsRhoJsNow
    int diag_flag;

    int nrequests;

    //! Tells which iteration was last time the patches moved (by moving window or load balancing)
    unsigned int lastIterationPatchesMoved;

    DomainDecomposition* domain_decomposition_;
        

    //! Methods to access readably to patch PIC operators.
    //!   - patches_ should not be access outsied of VectorPatch
    //!   - for now in SimWindow
    inline Species* species(int ipatch, int ispec) {
        return (*this)(ipatch)->vecSpecies[ispec];
    }

    inline ElectroMagn* emfields(int ipatch) {
        return (*this)(ipatch)->EMfields;
    }

    inline Interpolator* interp(int ipatch){
        return (*this)(ipatch)->Interp;
    }

    inline Projector* proj(int ipatch){
        return (*this)(ipatch)->Proj;
    }

    inline PartWalls* partwalls(int ipatch){
        return (*this)(ipatch)->partWalls;
    }

private :
    
    //  Internal balancing members
    // ---------------------------
    std::vector<Patch*> recv_patches_;

    std::vector<int> recv_patch_id_;
    std::vector<int> send_patch_id_;

    //! Current intensity of antennas
    double antenna_intensity;

    std::vector<Timer*> diag_timers;
    
};


#endif
