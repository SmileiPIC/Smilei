#ifndef DOMAIN_H
#define DOMAIN_H

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
    void clean();
    
    DomainDecomposition* cartGeom_;
    Patch* cartPatch_;
    VectorPatch VecPatchCart_;
    Diagnostic* diagCart_; 
   
};

#endif
