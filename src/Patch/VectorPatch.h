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
class SimWindow; 

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

 

    void runAllDiags(Params& params, SmileiMPI* smpi, int* diag_flag, int itime, std::vector<Timer>& timer);
    void dynamics(Params& params, SmileiMPI* smpi, SimWindow* simWindow, int* diag_flag, double time_dual, std::vector<Timer>& timer);
    void sumDensities( int* diag_flag, std::vector<Timer>& timer );
    void solveMaxwell(Params& params, SimWindow* simWindow, int itime, double time_dual, std::vector<Timer>& timer);




    void createPatches(Params& params, SmileiMPI* smpi, SimWindow* simWindow);
    void setNbrParticlesToExch(SmileiMPI* smpi);
    //void exchangePatches(SmileiMPI* smpi);
    void exchangePatches(SmileiMPI* smpi, Params& params);
    void output_exchanges(SmileiMPI* smpi);
    



    void solvePoisson( Params &params, SmileiMPI* smpi );
    bool isRhoNull( SmileiMPI* smpi );


    void clear() {patches_.clear();}

    std::vector<Patch*> patches_;
    std::vector<Patch*> recv_patches_;

    std::vector<int> recv_patch_id_;
    std::vector<int> send_patch_id_;

    Diagnostic* Diags;

 private :
    // 1st patch index of patches_ (stored for balancing op)
    int refHindex_;

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

    inline std::vector<PartWall*> partwalls(int ipatch){
	return (*this)(ipatch)->vecPartWall;
    }


};


#endif
