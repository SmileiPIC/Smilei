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
   
};

#endif
