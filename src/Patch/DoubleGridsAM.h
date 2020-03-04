
#ifndef DOUBLEGRIDSAM_H
#define DOUBLEGRIDSAM_H

class Region;
class VectorPatch;
class Patch;
class Params;
class SmileiMPI;
class Timers;
class Field;
class ElectroMagn;
class ElectroMagnAM;

class DoubleGridsAM
{
public :

    static void syncCurrentsOnRegion( VectorPatch &vecPatches, Region &region, Params &params, SmileiMPI *smpi, Timers &timers, int itime, unsigned int imode );
    static void currentsOnRegionSend( ElectroMagnAM* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI * smpi, Patch *patch, Params& params, unsigned int imode );
    static void currentsOnRegionSendFinalize( ElectroMagnAM* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI* smpi, Patch *patch, Params& params, unsigned int imode );
    static void currentsOnRegionRecv( ElectroMagnAM* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Region& region, unsigned int imode );


    static void syncFieldsOnPatches( Region &region, VectorPatch &vecPatches, Params &params, SmileiMPI *smpi, Timers &timers, int itime, unsigned int imode );
    static void fieldsOnPatchesRecv( ElectroMagnAM* localfields, unsigned int hindex, int recv_from_global_patch_rank, SmileiMPI* smpi, Patch* patch, Params& params, unsigned int imode );
    static void fieldsOnPatchesRecvFinalize( ElectroMagnAM* localfields, unsigned int hindex, int recv_from_global_patch_rank, SmileiMPI* smpi, Patch* patch, unsigned int imode );
    static void fieldsOnPatchesSend( ElectroMagnAM* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Region& region, unsigned int imode );

    static void syncBOnPatches( Region &region, VectorPatch &vecPatches, Params &params, SmileiMPI *smpi, Timers &timers, int itime, unsigned int imode );
    static void bOnPatchesRecv( ElectroMagnAM* localfields, unsigned int hindex, int recv_from_global_patch_rank, SmileiMPI* smpi, Patch* patch, Params& params, unsigned int imode );
    static void bOnPatchesRecvFinalize( ElectroMagnAM* localfields, unsigned int hindex, int recv_from_global_patch_rank, SmileiMPI* smpi, Patch* patch, unsigned int imode );
    static void bOnPatchesSend( ElectroMagnAM* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Region& region, unsigned int imode );


    static void syncFieldsOnRegion( VectorPatch& vecPatches, Region& region, Params &params, SmileiMPI* smpi, unsigned int imode );
    static void fieldsOnRegionSend( ElectroMagnAM* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI * smpi, Patch *patch, Params& params, unsigned int imode );
    static void fieldsOnRegionSendFinalize( ElectroMagnAM* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI* smpi, Patch *patch, Params& params, unsigned int imode );
    static void fieldsOnRegionRecv( ElectroMagnAM* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Region& region, unsigned int imode  );

    static void syncCurrentsOnPatches( Region &region, VectorPatch &vecPatches, Params &params, SmileiMPI *smpi, Timers &timers, int itime, unsigned int imode );
    static void currentsOnPatchesRecv( ElectroMagnAM* localfields, unsigned int hindex, int recv_from_global_patch_rank, SmileiMPI* smpi, Patch* patch, Params& params, unsigned int imode );
    static void currentsOnPatchesRecvFinalize( ElectroMagnAM* localfields, unsigned int hindex, int recv_from_global_patch_rank, SmileiMPI* smpi, Patch* patch, unsigned int imode );
    static void currentsOnPatchesSend( ElectroMagnAM* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Region& region, unsigned int imode );
};

#endif
