
#ifndef DOUBLEGRIDS_H
#define DOUBLEGRIDS_H

class Region;
class VectorPatch;
class Patch;
class Params;
class SmileiMPI;
class Timers;
class Field;
class ElectroMagn;

class DoubleGrids
{
public :

    static void syncCurrentsOnRegion( VectorPatch &vecPatches, Region &region, Params &params, SmileiMPI *smpi, Timers &timers );
    static void currentsOnRegionSend( ElectroMagn* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI * smpi, Patch *patch, Params& params );
    static void currentsOnRegionSendFinalize( Patch *patch, Params& params );
    static void currentsOnRegionRecv( ElectroMagn* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Region& region );

    static void syncFieldsOnPatches( Region &region, VectorPatch &vecPatches, Params &params, SmileiMPI *smpi, Timers &timers );
    static void fieldsOnPatchesRecv( ElectroMagn* localfields, unsigned int hindex, int recv_from_global_patch_rank, SmileiMPI* smpi, Patch* patch );
    static void fieldsOnPatchesRecvFinalize( Patch* patch );
    static void fieldsOnPatchesSend( ElectroMagn* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Region& region );

    static void syncBOnPatches( Region &region, VectorPatch &vecPatches, Params &params, SmileiMPI *smpi, Timers &timers );
    static void bOnPatchesRecv( ElectroMagn* localfields, unsigned int hindex, int recv_from_global_patch_rank, SmileiMPI* smpi, Patch* patch );
    static void bOnPatchesRecvFinalize( Patch* patch );
    static void bOnPatchesSend( ElectroMagn* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Region& region );

    static void syncFieldsOnRegion( VectorPatch& vecPatches, Region& region, Params &params, SmileiMPI* smpi );
    static void fieldsOnRegionSend( ElectroMagn* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI * smpi, Patch *patch, Params& params );
    static void fieldsOnRegionSendFinalize( Patch* patch );
    static void fieldsOnRegionRecv( ElectroMagn* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Region& region );

    static void syncCurrentsOnPatches( Region &region, VectorPatch &vecPatches, Params &params, SmileiMPI *smpi, Timers &timers );
    static void currentsOnPatchesRecv( ElectroMagn* localfields, unsigned int hindex, int recv_from_global_patch_rank, SmileiMPI* smpi, Patch* patch );
    static void currentsOnPatchesRecvFinalize( Patch* patch );
    static void currentsOnPatchesSend( ElectroMagn* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Region& region );

};

#endif
