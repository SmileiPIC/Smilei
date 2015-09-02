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

#include "DiagParams.h"
#include "PicParams.h"
#include "LaserParams.h"
#include "SmileiMPI.h"
#include "SimWindow.h"
#include "Diagnostic.h"
#include "SmileiIO.h"

class Diagnostic;
class DiagnosticScalar;
class Field;

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

    void exchangeParticles(int ispec, PicParams &params, SmileiMPI* smpi);
    void exchangeParticles(int ispec, PicParams &params);
    void sumRhoJ( unsigned int diag_flag );
    void sumRhoJs( int ispec );
    void exchangeE(  );
    void exchangeB(  );

    void computeGlobalDiags(int timestep);
    void computeScalarsDiags(int timestep);
    void computePhaseSpace();

    void initProbesDiags(PicParams& params, DiagParams &diag_params, int timestep);
    void finalizeProbesDiags(PicParams& params, DiagParams &diag_params, int timestep);
    void definePatchDiagsMaster();

    void createPatches(PicParams& params, DiagParams& diag_params, LaserParams& laser_params, SmileiMPI* smpi, SimWindow* simWindow);
    void setNbrParticlesToExch(SmileiMPI* smpi);
    void exchangePatches(SmileiMPI* smpi);
    void exchangePatches_new(SmileiMPI* smpi);
    
    void initDumpFields(PicParams& params, DiagParams &diag_params, int timestep);
    void finalizeDumpFields(PicParams& params, DiagParams &diag_params, int timestep);


    void solvePoisson( PicParams &params, SmileiMPI* smpi );
    bool isRhoNull( SmileiMPI* smpi );

    void exchange( std::vector<Field*> fields );
    void exchange0( std::vector<Field*> fields );
    void exchange1( std::vector<Field*> fields );
    void sum( std::vector<Field*> fields );


    void clear() {patches_.clear();}

    std::vector<Patch*> patches_;
    std::vector<Patch*> recv_patches_;

    std::vector<int> recv_patch_id_;
    std::vector<int> send_patch_id_;    

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
