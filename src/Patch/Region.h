#ifndef REGION_H
#define REGION_H
#include "interface.h"
#include "VectorPatch.h"

class DomainDecomposition;
class Patch;
class Diagnostic;

class Params;
class SimWindow;
class Timers;

class Region
{
public:
    Region( Params &params );
    ~Region();
    
    void build( Params &params, SmileiMPI *smpi, VectorPatch &vecPatches, OpenPMDparams &openPMD, bool global_region );
    void coupling( Params &params, bool global_region );
    void solveMaxwell( Params &params, SimWindow *simWindow, int itime, double time_dual, Timers &timers, SmileiMPI *smpi );
    void solveEnvelope( Params &params, SimWindow *simWindow, int itime, double time_dual, Timers &timers, SmileiMPI *smpi );
    void clean();
    VectorPatch vecPatch_;
    DomainDecomposition *decomposition_;
    Patch *patch_;
    Diagnostic *diag_;
    
    int hrank_global_region( int hindex, Params& params, VectorPatch& vecPatches );

    void identify_additional_patches(SmileiMPI* smpi, VectorPatch& vecPatches, Params& params, SimWindow* simWindow);
    std::vector<int> additional_patches_;
    std::vector<int> additional_patches_ranks;
    std::vector<int> local_patches_;
    void identify_missing_patches(SmileiMPI* smpi, VectorPatch& vecPatches, Params& params);
    std::vector<int> missing_patches_;
    std::vector<int> missing_patches_ranks;

    void reset_fitting(SmileiMPI* smpi, Params& params);

    void reset_mapping();

    Patch* fake_patch;
private:
    bool coupled_;
};

#endif
