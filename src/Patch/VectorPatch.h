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

#include "Params.h"
#include "LaserParams.h"
#include "SmileiMPI.h"
#include "SimWindow.h"
#include "Diagnostic.h"
#include "SmileiIO.h"

class Diagnostic;
class DiagnosticScalar;
class Field;
class Timer;

//! Class Patch : sub MPI domain 
//!     Collection of patch = MPI domain

class VectorPatch {
 public :
    VectorPatch();
    ~VectorPatch();

    void resize(int npatches) {patches_.resize(npatches);};
    int size() const {return patches_.size();};

    inline Patch* operator()(int ipatch) {return patches_[ipatch];};
    inline void set_refHindex() {refHindex_ = patches_[0]->Hindex();};

    void exchangeParticles(int ispec, Params &params, SmileiMPI* smpi);
#ifdef _NOTFORNOW
    void exchangeParticles(int ispec, Params &params);
#endif
    void sumRhoJ( unsigned int diag_flag );
    void sumRhoJs( int ispec );
    void exchangeE(  );
    void exchangeB(  );

    void runAllDiags(Params& params, SmileiMPI* smpi, int* diag_flag, int itime, std::vector<Timer>& timer);

    void computeGlobalDiags(int timestep);
    void computeScalarsDiags(int timestep);
    void computePhaseSpace();
    void computeParticlesDiags(int timestep);

    void initProbesDiags(Params& params, int timestep);
    void finalizeProbesDiags(Params& params, int timestep);
    void definePatchDiagsMaster(hid_t globalFile, hid_t globalFileAvg);
    void definePatchDiagsMaster();
    void updatePatchFieldDump( Params& params );

    void createPatches(Params& params, SmileiMPI* smpi, SimWindow* simWindow);
    void setNbrParticlesToExch(SmileiMPI* smpi);
    //void exchangePatches(SmileiMPI* smpi);
    void exchangePatches(SmileiMPI* smpi, Params& params);
    void output_exchanges(SmileiMPI* smpi);
    
    void initDumpFields(Params& params, int timestep);
    void finalizeDumpFields(Params& params, int timestep);

    void initTrackParticles(Params& params, SmileiMPI* smpi);

    void initCollisionDebug();


    void solvePoisson( Params &params, SmileiMPI* smpi );
    bool isRhoNull( SmileiMPI* smpi );

    void exchange( std::vector<Field*> fields );
    void exchange0( std::vector<Field*> fields );
    void exchange1( std::vector<Field*> fields );
    void sum( std::vector<Field*> fields );


    void clear() {patches_.clear();}
    void resizeFields();

    std::vector<Patch*> patches_;
    std::vector<Patch*> recv_patches_;

    std::vector<int> recv_patch_id_;
    std::vector<int> send_patch_id_;

    Diagnostic* Diags;

 private :
    // 1st patch index of patches_ (stored for balancing op)
    int refHindex_;

    std::vector<Field*> Ap_;
    std::vector<Field*> Bx_;
    std::vector<Field*> By_;
    std::vector<Field*> Bz_;
    std::vector<Field*> Ex_;
    std::vector<Field*> Ey_;
    std::vector<Field*> Ez_;

    std::vector<Field*> Jx_;
    std::vector<Field*> Jy_;
    std::vector<Field*> Jz_;
    std::vector<Field*> rho_;

};


#endif
