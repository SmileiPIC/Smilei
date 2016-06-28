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
#include "DiagnosticFactory.h"

#include "Params.h"
#include "SmileiMPI.h"
#include "SimWindow.h"

class Field;
class Timer;
class SimWindow; 

//! Class Patch : sub MPI domain 
//!     Collection of patch = MPI domain

class VectorPatch {
public :

    VectorPatch();
    ~VectorPatch();
    
    void close(SmileiMPI*);

    //! VectorPatch = 
    //! - std::vector<Patch*>
    //! - interfaces between main programs & main PIC operators
    //! - methods to balance computation
    std::vector<Patch*> patches_;

    std::vector<Diagnostic*> globalDiags;
    std::vector<Diagnostic*> localDiags;


    //! Some vector operations extended to VectorPatch
    inline void resize(int npatches) {
        patches_.resize(npatches);
    }
    inline int  size() const {
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

    void createDiags(Params& params, SmileiMPI* smpi);

    //! get a particular scalar
    inline double getScalar(std::string name) {
        DiagnosticScalar* diag = static_cast<DiagnosticScalar*>( globalDiags[0] );
        return diag->getScalar( name );
    }

    bool fieldTimeIsNow( int timestep ) {
        if( fieldsTimeSelection!=NULL )
            return fieldsTimeSelection->theTimeIsNow(timestep);
        else
            return false;
    }
    
    bool printScalars( int timestep ) {
        return static_cast<DiagnosticScalar*>(globalDiags[0])->printNow(timestep);
    }


   
    // Interfaces between main programs & main PIC operators
    // -----------------------------------------------------

    //! For all patch, move particles (restartRhoJ(s), dynamics and exchangeParticles)
    void dynamics(Params& params, SmileiMPI* smpi, SimWindow* simWindow, int* diag_flag, double time_dual,
                  std::vector<Timer>& timer);

    //! For all patch, sum densities on ghost cells (sum per species if needed, sync per patch and MPI sync)
    void sumDensities( int* diag_flag, std::vector<Timer>& timer );

    //! For all patch, update E and B (Ampere, Faraday, boundary conditions, exchange B and center B)
    void solveMaxwell(Params& params, SimWindow* simWindow, int itime, double time_dual,
                      std::vector<Timer>& timer);

    //! For all patch, Compute and Write all diags (Scalars, Probes, Phases, TrackParticles, Fields, Average fields)
    void runAllDiags(Params& params, SmileiMPI* smpi, int* diag_flag, int itime, std::vector<Timer>& timer);
    void initAllDiags(Params& params, SmileiMPI* smpi);
    void closeAllDiags(SmileiMPI* smpi);
    void openAllDiags(Params& params, SmileiMPI* smpi);

    //! Check if rho is null (MPI & patch sync)
    bool isRhoNull( SmileiMPI* smpi );

    //! Solve Poisson to initialize E
    void solvePoisson( Params &params, SmileiMPI* smpi );
    
    //! For all patch initialize the externals (lasers, fields, antennas)
    void initExternals(Params& params);



    //  Balancing methods
    // ------------------

    //! Explicits patch movement regarding new patch distribution stored in smpi->patch_count
    void createPatches(Params& params, SmileiMPI* smpi, SimWindow* simWindow);

    //! Exchange patches, based on createPatches initialization
    void exchangePatches(SmileiMPI* smpi, Params& params);

    //! Write in a file patches communications
    void output_exchanges(SmileiMPI* smpi);

    // Lists of fields
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
    bool hasAntennas;
    
    //! 1st patch index of patches_ (stored for balancing op)
    int refHindex_;
    
    //! Copy of the fields time selection
    TimeSelection * fieldsTimeSelection;

 private :

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

    //  Internal balancing members
    // ---------------------------
    std::vector<Patch*> recv_patches_;

    std::vector<int> recv_patch_id_;
    std::vector<int> send_patch_id_;

    
};


#endif
