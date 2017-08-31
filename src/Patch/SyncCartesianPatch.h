
#ifndef SYNCCARTESIANPATCH_H
#define SYNCCARTESIANPATCH_H

class VectorPatch;
class Patch;
class Params;
class SmileiMPI;
class Timers;
class Field;

class SyncCartesianPatch {
public :

    static void patchedToCartesian( VectorPatch& vecPatches, Patch* patch, Params &params, SmileiMPI* smpi, Timers &timers, int itime );
    static void sync( Field* inField, Field* outField, Params &params, SmileiMPI* smpi );

};

#endif
