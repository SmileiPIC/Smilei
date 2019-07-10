
#ifndef DOUBLEGRIDS_H
#define DOUBLEGRIDS_H

class Domain;
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

    static void syncCurrentsOnDomain( VectorPatch &vecPatches, Domain &domain, Params &params, SmileiMPI *smpi, Timers &timers, int itime );
    static void currentsOnDomainSend( ElectroMagn* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI * smpi, Patch *patch, Params& params );
    static void currentsOnDomainSendFinalize( ElectroMagn* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI* smpi, Patch *patch, Params& params );
    static void currentsOnDomainRecv( ElectroMagn* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Domain& domain );

    static void syncFieldsOnPatches( Domain &domain, VectorPatch &vecPatches, Params &params, SmileiMPI *smpi, Timers &timers, int itime );
    static void fieldsOnPatchesRecv( ElectroMagn* localfields, unsigned int hindex, int recv_from_global_patch_rank, SmileiMPI* smpi, Patch* patch );
    static void fieldsOnPatchesRecvFinalize( ElectroMagn* localfields, unsigned int hindex, int recv_from_global_patch_rank, SmileiMPI* smpi, Patch* patch );
    static void fieldsOnPatchesSend( ElectroMagn* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Domain& domain );

    static void syncFieldsOnDomain( VectorPatch& vecPatches, Domain& domain, Params &params, SmileiMPI* smpi );
    static void fieldsOnDomainSend( ElectroMagn* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI * smpi, Patch *patch, Params& params );
    static void fieldsOnDomainSendFinalize( ElectroMagn* localfields, unsigned int hindex, int send_to_global_patch_rank, SmileiMPI* smpi, Patch *patch, Params& params );
    static void fieldsOnDomainRecv( ElectroMagn* globalfields, unsigned int hindex, int local_patch_rank, VectorPatch& vecPatches, Params &params, SmileiMPI* smpi, Domain& domain );

};

#endif
