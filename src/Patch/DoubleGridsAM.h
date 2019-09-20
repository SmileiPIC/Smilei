
#ifndef DOUBLEGRIDSAM_H
#define DOUBLEGRIDSAM_H

class Domain;
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

    static void syncCurrentsOnDomain( VectorPatch &vecPatches, Domain &domain, Params &params, SmileiMPI *smpi, Timers &timers, int itime, unsigned int imode );
    static void currentsOnDomainSend( ElectroMagnAM* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI * smpi, Patch *patch, Params& params, unsigned int imode );
    static void currentsOnDomainSendFinalize( ElectroMagnAM* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI* smpi, Patch *patch, Params& params, unsigned int imode );
    static void currentsOnDomainRecv( ElectroMagnAM* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Domain& domain, unsigned int imode );


    static void syncFieldsOnPatches( Domain &domain, VectorPatch &vecPatches, Params &params, SmileiMPI *smpi, Timers &timers, int itime, unsigned int imode );
    static void fieldsOnPatchesRecv( ElectroMagnAM* localfields, unsigned int hindex, int recv_from_global_patch_rank, SmileiMPI* smpi, Patch* patch, Params& params, unsigned int imode );
    static void fieldsOnPatchesRecvFinalize( ElectroMagnAM* localfields, unsigned int hindex, int recv_from_global_patch_rank, SmileiMPI* smpi, Patch* patch, unsigned int imode );
    static void fieldsOnPatchesSend( ElectroMagnAM* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Domain& domain, unsigned int imode );

    static void syncBOnPatches( Domain &domain, VectorPatch &vecPatches, Params &params, SmileiMPI *smpi, Timers &timers, int itime, unsigned int imode );
    static void bOnPatchesRecv( ElectroMagnAM* localfields, unsigned int hindex, int recv_from_global_patch_rank, SmileiMPI* smpi, Patch* patch, Params& params, unsigned int imode );
    static void bOnPatchesRecvFinalize( ElectroMagnAM* localfields, unsigned int hindex, int recv_from_global_patch_rank, SmileiMPI* smpi, Patch* patch, unsigned int imode );
    static void bOnPatchesSend( ElectroMagnAM* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Domain& domain, unsigned int imode );


    static void syncFieldsOnDomain( VectorPatch& vecPatches, Domain& domain, Params &params, SmileiMPI* smpi, unsigned int imode );
    static void fieldsOnDomainSend( ElectroMagnAM* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI * smpi, Patch *patch, Params& params, unsigned int imode );
    static void fieldsOnDomainSendFinalize( ElectroMagnAM* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI* smpi, Patch *patch, Params& params, unsigned int imode );
    static void fieldsOnDomainRecv( ElectroMagnAM* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Domain& domain, unsigned int imode  );

};

#endif
