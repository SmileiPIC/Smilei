#ifndef DOMAIN_H
#define DOMAIN_H
#include "interface.h"
#include "VectorPatch.h" 

class DomainDecomposition;
class Patch;
class Diagnostic;

class Params;
class SimWindow;
class Timers;

class Domain
{
public:
    Domain( Params& params );
    ~Domain();

    void build( Params& params, SmileiMPI* smpi, VectorPatch& vecPatches, OpenPMDparams& openPMD );
    void solveMaxwell( Params& params, SimWindow* simWindow, int itime, double time_dual, Timers& timers );
    void solveEnvelope( Params& params, SimWindow* simWindow, int itime, double time_dual, Timers& timers );
    void clean();
    VectorPatch vecPatch_;
    DomainDecomposition* decomposition_;
    Patch* patch_;
    Diagnostic* diag_; 
   
    int hrank_global_domain( int hindex, Params& params, VectorPatch& vecPatches );

    void identify_additional_patches(SmileiMPI* smpi, VectorPatch& vecPatches, Params& params, SimWindow* simWindow);
    std::vector<int> additional_patches_;
    std::vector<int> additional_patches_ranks;
    std::vector<int> local_patches_;
    void identify_missing_patches(SmileiMPI* smpi, VectorPatch& vecPatches, Params& params);
    std::vector<int> missing_patches_;
    std::vector<int> missing_patches_ranks;

    void reset_mapping();

    Patch* fake_patch;



};

#endif
