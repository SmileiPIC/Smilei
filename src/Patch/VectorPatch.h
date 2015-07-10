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

//! Class Patch : sub MPI domain 
//!     Collection of patch = MPI domain

class VectorPatch {
 public :
    VectorPatch();
    ~VectorPatch();

    void resize(int npatches) {patches_.resize(npatches);};
    int size() const {return patches_.size();};

    inline Patch* operator()(int ipatch) {return patches_[ipatch];};

    void exchangeParticles(int ispec, PicParams &params, SmileiMPI* smpi);
    void sumRhoJ( unsigned int diag_flag );
    void sumRhoJs( int ispec );
    void exchangeE(  );
    void exchangeB(  );

    void computeGlobalDiags(int timestep);
    void computeScalarsDiags(int timestep);

    void initProbesDiags(PicParams& params, DiagParams &diag_params, int timestep);
    void finalizeProbesDiags(PicParams& params, DiagParams &diag_params, int timestep);
    void definePatchDiagsMaster();

    void createPatches(PicParams& params, DiagParams& diag_params, LaserParams& laser_params, SmileiMPI* smpi, SimWindow* simWindow);
    void setNbrParticlesToExch(SmileiMPI* smpi);
    void exchangePatches(SmileiMPI* smpi);
    
    void initDumpFields(PicParams& params, DiagParams &diag_params, int timestep);
    void finalizeDumpFields(PicParams& params, DiagParams &diag_params, int timestep);


    void solvePoisson( PicParams &params, SmileiMPI* smpi );


    void clear() {patches_.clear();}

    std::vector<Patch*> patches_;
    std::vector<Patch*> recv_patches_;

    std::vector<int> recv_patch_id_;
    std::vector<int> send_patch_id_;    

 private :
    // 1st patch index of patches_ (stored for balancing op)
    int refHindex_;
};


#endif
