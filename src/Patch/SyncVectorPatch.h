
#ifndef SYNCVECTORPATCH_H
#define SYNCVECTORPATCH_H

#include <vector>

class VectorPatch;
class Params;
class SmileiMPI;
class Field;
class Timers;

class SyncVectorPatch {
public :

    //! Particles synchronization
    static void exchangeParticles(VectorPatch& vecPatches, int ispec, Params &params, SmileiMPI* smpi, Timers &timers, int itime);
    static void finalize_and_sort_parts(VectorPatch& vecPatches, int ispec, Params &params, SmileiMPI* smpi, Timers &timers, int itime);

    //! Densities synchronization
    static void sumRhoJ  ( Params& params, VectorPatch& vecPatches, SmileiMPI* smpi, Timers &timers, int itime );
    static void sumRhoJs ( Params& params, VectorPatch& vecPatches, int ispec, SmileiMPI* smpi, Timers &timers, int itime );
    static void sum      ( std::vector<Field*> fields, VectorPatch& vecPatches, SmileiMPI* smpi, Timers &timers, int itime );
    static void sum_all_components ( std::vector<Field*>& fields, VectorPatch& vecPatches, SmileiMPI* smpi, Timers &timers, int itime );

    //! Fields synchronization
    static void exchangeE( Params& params, VectorPatch& vecPatches, SmileiMPI* smpi );
    static void finalizeexchangeE( Params& params, VectorPatch& vecPatches );
    static void exchangeB( Params& params, VectorPatch& vecPatches, SmileiMPI* smpi );
    static void finalizeexchangeB( Params& params, VectorPatch& vecPatches );
    static void exchangeJ( Params& params, VectorPatch& vecPatches, SmileiMPI* smpi );
    static void finalizeexchangeJ( Params& params, VectorPatch& vecPatches );

    static void exchange_along_all_directions         ( std::vector<Field*> fields, VectorPatch& vecPatches, SmileiMPI* smpi );
    static void finalize_exchange_along_all_directions( std::vector<Field*> fields, VectorPatch& vecPatches );
    static void exchange_along_all_directions_noomp         ( std::vector<Field*> fields, VectorPatch& vecPatches, SmileiMPI* smpi  );
    static void finalize_exchange_along_all_directions_noomp( std::vector<Field*> fields, VectorPatch& vecPatches );
    static void exchange_synchronized_per_direction   ( std::vector<Field*> fields, VectorPatch& vecPatches, SmileiMPI* smpi );

    static void exchange_all_components_along_X          ( std::vector<Field*>& fields, VectorPatch& vecPatches, SmileiMPI* smpi );
    static void finalize_exchange_all_components_along_X ( std::vector<Field*>& fields, VectorPatch& vecPatches );
    static void exchange_all_components_along_Y          ( std::vector<Field*>& fields, VectorPatch& vecPatches, SmileiMPI* smpi );
    static void finalize_exchange_all_components_along_Y ( std::vector<Field*>& fields, VectorPatch& vecPatches );
    static void exchange_all_components_along_Z          ( std::vector<Field*> fields, VectorPatch& vecPatches, SmileiMPI* smpi );
    static void finalize_exchange_all_components_along_Z ( std::vector<Field*> fields, VectorPatch& vecPatches );

    //! Deprecated field functions
    static void exchange_along_X          ( std::vector<Field*> fields, VectorPatch& vecPatches, SmileiMPI* smpi );
    static void finalize_exchange_along_X ( std::vector<Field*> fields, VectorPatch& vecPatches );
    static void exchange_along_Y          ( std::vector<Field*> fields, VectorPatch& vecPatches, SmileiMPI* smpi );
    static void finalize_exchange_along_Y ( std::vector<Field*> fields, VectorPatch& vecPatches );
    static void exchange_along_Z          ( std::vector<Field*> fields, VectorPatch& vecPatches, SmileiMPI* smpi );
    static void finalize_exchange_along_Z ( std::vector<Field*> fields, VectorPatch& vecPatches );

};

#endif
