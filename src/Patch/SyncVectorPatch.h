
#ifndef SYNCVECTORPATCH_H
#define SYNCVECTORPATCH_H

#include <vector>

class VectorPatch;
class Params;
class SmileiMPI;
class Field;

class SyncVectorPatch {
public :

    static void exchangeParticles(VectorPatch& vecPatches, int ispec, Params &params, SmileiMPI* smpi);
    static void sumRhoJ  ( VectorPatch& vecPatches, unsigned int diag_flag );
    static void sumRhoJs ( VectorPatch& vecPatches, int ispec );
    static void exchangeE( VectorPatch& vecPatches );
    static void exchangeB( VectorPatch& vecPatches );
    static void sum      ( std::vector<Field*> fields, VectorPatch& vecPatches );
    static void exchange ( std::vector<Field*> fields, VectorPatch& vecPatches );
    static void exchange0( std::vector<Field*> fields, VectorPatch& vecPatches );
    static void exchange1( std::vector<Field*> fields, VectorPatch& vecPatches );

};

#endif
