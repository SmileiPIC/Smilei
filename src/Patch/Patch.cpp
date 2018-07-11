
// CLOBAL COORDINATES:
//                           Patch_minGlobal                                                                      Patch_maxGlobal
//                      --------<===================================== gs ===================================>------------
//     GLOBAL INDICES:          0                                  .                                        nspace_global
//                           ix+oversize                                                                  ix+oversize
//                      ------------------------------------       .              ------------------------------------
//                      |   |   |     ...          |   |   |       .              |   |   |   |   ...    |   |   |   |
//                      |   |   |     ...          |   |   |       .              |   |   |   |   ...    |   |   |   |
//                      ------------------------------------       .              ------------------------------------
//                          Patch_minLocal    Patch_maxLocal       .             Patch_minLocal        Patch_maxLocal
//                                                 ----------------------------------------
//                                                 |   |   |       .              |   |   |
//                                                 |   |   |       .              |   |   |
//                                                 ----------------------------------------
// LOCAL COORDINATES:                             x(0) rlb        x(ix)             rub  x(nspace)
//                                                 ----<============= length =========>----
//     LOCAL INDICES:                              0   lb                            ub   nspace

#include "Patch.h"

#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "Hilbert_functions.h"
#include "PatchesFactory.h"
#include "SpeciesFactory.h"
#include "Particles.h"
#include "ElectroMagnFactory.h"
#include "ElectroMagnBC_Factory.h"
#include "DiagnosticFactory.h"
#include "CollisionsFactory.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// Patch constructor :
//   Called by PatchXD constructor which will finalize initialization
// ---------------------------------------------------------------------------------------------------------------------
Patch::Patch(Params& params, SmileiMPI* smpi, DomainDecomposition* domain_decomposition, unsigned int ipatch, unsigned int n_moved)
{

    hindex = ipatch;
    nDim_fields_ = params.nDim_field;

    initStep1(params);

#ifdef  __DETAILED_TIMERS
    // Initialize timers
    patch_timers.resize(13,0.);
#endif

} // END Patch::Patch




// Cloning patch constructor
Patch::Patch(Patch* patch, Params& params, SmileiMPI* smpi, DomainDecomposition* domain_decomposition, unsigned int ipatch, unsigned int n_moved, bool with_particles = true) {

    hindex = ipatch;
    nDim_fields_ = patch->nDim_fields_;

    initStep1(params);

#ifdef  __DETAILED_TIMERS
    // Initialize timers
    patch_timers.resize(13,0.);
#endif

}

void Patch::initStep1(Params& params)
{
    // for nDim_fields = 1 : bug if Pcoordinates.size = 1 !!
    //Pcoordinates.resize(nDim_fields_);
    Pcoordinates.resize( 2 );

    nbNeighbors_ = 2;
    neighbor_.resize(nDim_fields_);
    tmp_neighbor_.resize(nDim_fields_);
    send_tags_.resize(nDim_fields_);
    recv_tags_.resize(nDim_fields_);
    for ( int iDim = 0 ; iDim < nDim_fields_ ; iDim++ ) {
        neighbor_[iDim].resize(2,MPI_PROC_NULL);
        tmp_neighbor_[iDim].resize(2,MPI_PROC_NULL);
        send_tags_[iDim].resize(2,MPI_PROC_NULL);
        recv_tags_[iDim].resize(2,MPI_PROC_NULL);
    }
    MPI_neighbor_.resize(nDim_fields_);
    tmp_MPI_neighbor_.resize(nDim_fields_);
    for ( int iDim = 0 ; iDim < nDim_fields_; iDim++ ) {
        MPI_neighbor_[iDim].resize(2,MPI_PROC_NULL);
        tmp_MPI_neighbor_[iDim].resize(2,MPI_PROC_NULL);
    }

    oversize.resize( nDim_fields_ );
    for ( int iDim = 0 ; iDim < nDim_fields_; iDim++ )
        oversize[iDim] = params.oversize[iDim];
}


void Patch::initStep3( Params& params, SmileiMPI* smpi, unsigned int n_moved ) {
    // Compute MPI neighborood
    updateMPIenv(smpi);

    // Compute patch boundaries
    min_local.resize(params.nDim_field, 0.);
    max_local.resize(params.nDim_field, 0.);
    center   .resize(params.nDim_field, 0.);
    cell_starting_global_index.resize(params.nDim_field, 0);
    radius = 0.;
    for (unsigned int i = 0 ; i<params.nDim_field ; i++) {
        min_local[i] = ( Pcoordinates[i]   )*(params.n_space[i]*params.cell_length[i]);
        max_local[i] = ((Pcoordinates[i]+1))*(params.n_space[i]*params.cell_length[i]);
        cell_starting_global_index[i] += Pcoordinates[i]*params.n_space[i];
        cell_starting_global_index[i] -= params.oversize[i];
        center[i] = (min_local[i]+max_local[i])*0.5;
        radius += pow(max_local[i] - center[i] + params.cell_length[i], 2);
    }
    radius = sqrt(radius);

    cell_starting_global_index[0] += n_moved;
    min_local[0] += n_moved*params.cell_length[0];
    max_local[0] += n_moved*params.cell_length[0];
    center   [0] += n_moved*params.cell_length[0];

}

void Patch::finishCreation( Params& params, SmileiMPI* smpi, DomainDecomposition* domain_decomposition ) {
    // initialize vector of Species (virtual)
    vecSpecies = SpeciesFactory::createVector(params, this);

    // initialize the electromagnetic fields (virtual)
    EMfields   = ElectroMagnFactory::create(params, domain_decomposition, vecSpecies, this);

    // Create ad hoc interpolators and projectors for envelope
    if (params.Laser_Envelope_model){
        Interp_envelope      = InterpolatorFactory::create_env_interpolator(params, this, params.vectorization_mode == "normal");
    } // + patchId -> idx_domain_begin (now = ref smpi)
    else {
        Interp_envelope      = NULL;
    }

    // Initialize the collisions
    vecCollisions = CollisionsFactory::create(params, this, vecSpecies);

    // Initialize the particle walls
    partWalls = new PartWalls(params, this);

    // Initialize the probes
    probes = DiagnosticFactory::createProbes();

    probesInterp = InterpolatorFactory::create(params, this, false);

    if (has_an_MPI_neighbor())
        createType(params);

}


void Patch::finishCloning( Patch* patch, Params& params, SmileiMPI* smpi, unsigned int n_moved, bool with_particles = true ) {
    // clone vector of Species (virtual)
    vecSpecies = SpeciesFactory::cloneVector(patch->vecSpecies, params, this, with_particles);

    // clone the electromagnetic fields (virtual)
    EMfields   = ElectroMagnFactory::clone(patch->EMfields, params, vecSpecies, this, n_moved);

    // Create ad hoc interpolators and projectors for envelope
    if (params.Laser_Envelope_model){
        Interp_envelope  = InterpolatorFactory::create_env_interpolator(params, this, params.vectorization_mode == "normal");
    }

    // clone the collisions
    vecCollisions = CollisionsFactory::clone(patch->vecCollisions, params);

    // clone the particle walls
    partWalls = new PartWalls(patch->partWalls, this);

    // clone the probes
    probes = DiagnosticFactory::cloneProbes(patch->probes);

    probesInterp = InterpolatorFactory::create(params, this, false);

    if (has_an_MPI_neighbor())
        createType(params);

}

void Patch::finalizeMPIenvironment(Params& params) {
    int nb_comms(9); // E, B, B_m : min number of comms

    if (params.geometry == "3drz")
        nb_comms += 9*(params.nmodes - 1);
    // if envelope is present,
    // add to comms A, A0, Phi, Phi_old, GradPhi (x,y,z components), GradPhi_old (x,y,z components)
    if (params.Laser_Envelope_model){
        nb_comms += 10;
    }

    // add comms for species
    nb_comms += 2*vecSpecies.size();

    // Dynamic vectorization:
    if (params.has_dynamic_vectorization)
    {
        nb_comms += 2;
    }

    // Radiated energy
    if (params.hasMCRadiation ||
        params.hasLLRadiation ||
        params.hasNielRadiation)
    {
        for (int ispec=0 ; ispec<(int)this->vecSpecies.size() ; ispec++)
        {
            if (this->vecSpecies[ispec]->Radiate)
            {
                nb_comms ++;
            }
        }
    }

    // Just apply on species & fields to start

    for( unsigned int idiag=0; idiag<EMfields->allFields_avg.size(); idiag++) {
        nb_comms += EMfields->allFields_avg[idiag].size();
    }
    nb_comms += EMfields->antennas.size();

    for (unsigned int bcId=0 ; bcId<EMfields->emBoundCond.size() ; bcId++ ) {
        if(EMfields->emBoundCond[bcId]) {
            for (unsigned int laserId=0 ; laserId < EMfields->emBoundCond[bcId]->vecLaser.size() ; laserId++ ) {
                nb_comms += 4;
            }
        }
        if ( EMfields->extFields.size()>0 ) {
            if (dynamic_cast<ElectroMagnBC1D_SM*>(EMfields->emBoundCond[bcId]) )
                nb_comms += 4;
            else if ( dynamic_cast<ElectroMagnBC2D_SM*>(EMfields->emBoundCond[bcId]) )
                nb_comms += 12;
            else if ( dynamic_cast<ElectroMagnBC3D_SM*>(EMfields->emBoundCond[bcId]) )
                nb_comms += 18;
        }
    }
    requests_.resize( nb_comms, MPI_REQUEST_NULL );

}


void Patch::set( Params& params, DomainDecomposition* domain_decomposition, VectorPatch& vecPatch )
{
    Pcoordinates.resize( params.nDim_field );


    min_local = vecPatch(0)->min_local;
    max_local = vecPatch(0)->max_local;
    center   .resize( nDim_fields_, 0. );
    cell_starting_global_index = vecPatch(0)->cell_starting_global_index;
    radius = 0.;

    // Constraint to enforce 1 neighboor per side
    double nppp_root = pow( vecPatch.size(), 1./(double)nDim_fields_ );
    if ( fabs( (double)(int)nppp_root - nppp_root ) > 0. )
        ERROR( "Bad choice of decomposition" );

    for (int i = 0 ; i<nDim_fields_ ; i++) {

        for ( unsigned int ipatch = 0 ; ipatch < vecPatch.size() ; ipatch++  ) {
            if ( vecPatch(ipatch)->min_local[i] <= min_local[i] ) {
                min_local[i] = vecPatch(ipatch)->min_local[i];
                if (vecPatch(ipatch)->MPI_neighbor_[i][0]!=MPI_PROC_NULL)
                    MPI_neighbor_[i][0] = vecPatch(ipatch)->MPI_neighbor_[i][0];
                if (vecPatch(ipatch)->neighbor_[i][0]!=MPI_PROC_NULL)
                    neighbor_[i][0] = (vecPatch(ipatch)->neighbor_[i][0] / vecPatch.size() );
            }
            if ( vecPatch(ipatch)->max_local[i] >= max_local[i] ) {
                max_local[i] = vecPatch(ipatch)->max_local[i];
                if (vecPatch(ipatch)->MPI_neighbor_[i][1]!=MPI_PROC_NULL)
                    MPI_neighbor_[i][1] = vecPatch(ipatch)->MPI_neighbor_[i][1];
                if (vecPatch(ipatch)->neighbor_[i][1]!=MPI_PROC_NULL)
                    neighbor_[i][1] = (vecPatch(ipatch)->neighbor_[i][1] / vecPatch.size() );
            }
            if ( vecPatch(ipatch)->cell_starting_global_index[i] <= cell_starting_global_index[i] )
                cell_starting_global_index[i] = vecPatch(ipatch)->cell_starting_global_index[i];
        }
        Pcoordinates[i] = (cell_starting_global_index[i]+params.oversize[i]) / params.n_space[i] / params.global_factor[i];

        center[i] = (min_local[i]+max_local[i])*0.5;
        radius += pow(max_local[i] - center[i] + params.cell_length[i], 2);
    }
    radius = sqrt(radius);

    //cart_updateMPIenv(smpi);

    MPI_me_ = vecPatch(0)->MPI_me_;
    for ( int i = 0 ; i<nDim_fields_ ; i++) {
        if ((MPI_neighbor_[i][0]==MPI_me_) && (params.EM_BCs[i][0]!="periodic"))
            MPI_neighbor_[i][0] = MPI_PROC_NULL;
        if ((MPI_neighbor_[i][1]==MPI_me_) && (params.EM_BCs[i][0]!="periodic"))
            MPI_neighbor_[i][1] = MPI_PROC_NULL;
        if ((neighbor_[i][0]==MPI_me_) && (params.EM_BCs[i][0]!="periodic"))
            neighbor_[i][0] = MPI_PROC_NULL;
        if ((neighbor_[i][1]==MPI_me_) && (params.EM_BCs[i][0]!="periodic"))
            neighbor_[i][1] = MPI_PROC_NULL;
    }


    //cout << "MPI Nei\t"  << "\t" << MPI_neighbor_[1][1] << endl;
    //cout << "MPI Nei\t"  << MPI_neighbor_[0][0] << "\t" << MPI_me_ << "\t" << MPI_neighbor_[0][1] << endl;
    //cout << "MPI Nei\t"  << "\t" << MPI_neighbor_[1][0] << endl;


    for (int iDim=0 ; iDim< (int)neighbor_.size() ; iDim++)
        for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
            send_tags_[iDim][iNeighbor] = buildtag( hindex, iDim, iNeighbor, 5 );
            recv_tags_[iDim][iNeighbor] = buildtag( neighbor_[iDim][(iNeighbor+1)%2], iDim, iNeighbor, 5 );
        }

    EMfields   = ElectroMagnFactory::create(params, domain_decomposition, vecPatch(0)->vecSpecies, this);

    vecSpecies.resize(0);
    vecCollisions.resize(0);
    partWalls = NULL;
    probes.resize(0);

    if (has_an_MPI_neighbor())
        createType2(params);


}


// ---------------------------------------------------------------------------------------------------------------------
// Delete Patch members
// ---------------------------------------------------------------------------------------------------------------------
Patch::~Patch() {

    delete probesInterp;

    for(unsigned int i=0; i<probes.size(); i++)
        delete probes[i];

    for(unsigned int i=0; i<vecCollisions.size(); i++) delete vecCollisions[i];
    vecCollisions.clear();

    if (Interp_envelope     != NULL) delete Interp_envelope;

    if (partWalls!=NULL) delete partWalls;

    if (EMfields !=NULL) delete EMfields;

    for (unsigned int ispec=0 ; ispec<vecSpecies.size(); ispec++) delete vecSpecies[ispec];
    vecSpecies.clear();

} // END Patch::~Patch


// ---------------------------------------------------------------------------------------------------------------------
// Compute MPI rank of patch neigbors and current patch
// ---------------------------------------------------------------------------------------------------------------------
void Patch::updateTagenv(SmileiMPI* smpi)
{
    for (int iDim=0 ; iDim< (int)neighbor_.size() ; iDim++)
        for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
            send_tags_[iDim][iNeighbor] = buildtag( hindex, iDim, iNeighbor, 5 );
            recv_tags_[iDim][iNeighbor] = buildtag( neighbor_[iDim][(iNeighbor+1)%2], iDim, iNeighbor, 5 );
        }
}
void Patch::updateMPIenv(SmileiMPI* smpi)
{
    // HARDCODED VALUE - to test cartesian decomposition
//    vector<int> npattchpp(2,1); // NDOMAIN PER PROC - WORKS WITH 1 domain per process, ATTENTION PERIDODICITY
//    vector<int> npatchs(2,4);   // NDOMAIN TOTAL
//    vector<int> gCoord(2,0);
    MPI_me_ = smpi->smilei_rk;

    for (int iDim = 0 ; iDim < nDim_fields_ ; iDim++)
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++){
            MPI_neighbor_[iDim][iNeighbor] = smpi->hrank(neighbor_[iDim][iNeighbor]);

//            gCoord[0] = ( (int)Pcoordinates[0] + (1-iDim) * ((1-iNeighbor)*(-1) + (iNeighbor)) );
//            gCoord[1] = ( (int)Pcoordinates[1] + (iDim  ) * ((1-iNeighbor)*(-1) + (iNeighbor)) );
//            if ( gCoord[0] < 0 )
//                MPI_neighbor_[iDim][iNeighbor] = MPI_PROC_NULL;
//            else if ( gCoord[0] >= npatchs[0] )
//                MPI_neighbor_[iDim][iNeighbor] = MPI_PROC_NULL;
//            else if ( gCoord[1] < 0 )
//                MPI_neighbor_[iDim][iNeighbor] = MPI_PROC_NULL;
//            else if ( gCoord[1] >= npatchs[1] )
//                MPI_neighbor_[iDim][iNeighbor] = MPI_PROC_NULL;
//            else {
//                gCoord[0] = (double)( (int)Pcoordinates[0] + (1-iDim) * ((1-iNeighbor)*(-1) + (iNeighbor)) ) /(double)npattchpp[0];
//                gCoord[1] = (double)( (int)Pcoordinates[1] + (iDim  ) * ((1-iNeighbor)*(-1) + (iNeighbor)) ) /(double)npattchpp[1];
//                MPI_neighbor_[iDim][iNeighbor] = gCoord[0] * npatchs[1] + gCoord[1];
//            }
        }

    for (int iDim=0 ; iDim< (int)neighbor_.size() ; iDim++)
        for (int iNeighbor=0 ; iNeighbor<2 ; iNeighbor++) {
            send_tags_[iDim][iNeighbor] = buildtag( hindex, iDim, iNeighbor, 5 );
            recv_tags_[iDim][iNeighbor] = buildtag( neighbor_[iDim][(iNeighbor+1)%2], iDim, iNeighbor, 5 );
        }

} // END updateMPIenv

// ---------------------------------------------------------------------------------------------------------------------
// Clean the MPI buffers for communications
// ---------------------------------------------------------------------------------------------------------------------
void Patch::cleanMPIBuffers(int ispec, Params& params)
{
    int ndim = params.nDim_field;

    for (int iDim=0 ; iDim < ndim ; iDim++){
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
            vecSpecies[ispec]->MPIbuff.partRecv[iDim][iNeighbor].clear();//resize(0,ndim);
            vecSpecies[ispec]->MPIbuff.partSend[iDim][iNeighbor].clear();//resize(0,ndim);
            vecSpecies[ispec]->MPIbuff.part_index_send[iDim][iNeighbor].clear();
            //vecSpecies[ispec]->MPIbuff.part_index_send[iDim][iNeighbor].resize(0);
            vecSpecies[ispec]->MPIbuff.part_index_recv_sz[iDim][iNeighbor] = 0;
        }
    }
    vecSpecies[ispec]->indexes_of_particles_to_exchange.clear();
} // cleanMPIBuffers


// ---------------------------------------------------------------------------------------------------------------------
// Split particles Id to send in per direction and per patch neighbor dedicated buffers
// Apply periodicity if necessary
// ---------------------------------------------------------------------------------------------------------------------
void Patch::initExchParticles(SmileiMPI* smpi, int ispec, Params& params)
{
    Particles &cuParticles = (*vecSpecies[ispec]->particles);
    int ndim = params.nDim_field;
    int idim,check;
    std::vector<int>* indexes_of_particles_to_exchange = &vecSpecies[ispec]->indexes_of_particles_to_exchange;
//    double xmax[3];

    for (int iDim=0 ; iDim < ndim ; iDim++){
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
            vecSpecies[ispec]->MPIbuff.partRecv[iDim][iNeighbor].clear();//resize(0,ndim);
            vecSpecies[ispec]->MPIbuff.partSend[iDim][iNeighbor].clear();//resize(0,ndim);
            vecSpecies[ispec]->MPIbuff.part_index_send[iDim][iNeighbor].resize(0);
            vecSpecies[ispec]->MPIbuff.part_index_recv_sz[iDim][iNeighbor] = 0;
        }
    }

    int n_part_send = (*indexes_of_particles_to_exchange).size();

    int iPart;

    // Define where particles are going
    //Put particles in the send buffer it belongs to. Priority to lower dimensions.
    if ( params.geometry != "3drz") {
        for (int i=0 ; i<n_part_send ; i++) {
            iPart = (*indexes_of_particles_to_exchange)[i];
            check = 0;
            idim = 0;
            //Put indexes of particles in the first direction they will be exchanged and correct their position according to periodicity for the first exchange only.
            while (check == 0 && idim<ndim){
                if ( cuParticles.position(idim,iPart) < min_local[idim]){
                    if ( neighbor_[idim][0]!=MPI_PROC_NULL) {
                        vecSpecies[ispec]->MPIbuff.part_index_send[idim][0].push_back( iPart );
                    }
                    //If particle is outside of the global domain (has no neighbor), it will not be put in a send buffer and will simply be deleted.
                    check = 1;
                }
                else if ( cuParticles.position(idim,iPart) >= max_local[idim]){
                    if( neighbor_[idim][1]!=MPI_PROC_NULL) {
                        vecSpecies[ispec]->MPIbuff.part_index_send[idim][1].push_back( iPart );
                    }
                    check = 1;
                }
                idim++;
            }
        }
    }
    else { //if (geometry == "3drz")
        for (int i=0 ; i<n_part_send ; i++) {
            iPart = (*indexes_of_particles_to_exchange)[i];
            check = 0;
            idim = 0;
            //Put indexes of particles in the first direction they will be exchanged and correct their position according to periodicity for the first exchange only.
            while (check == 0 && idim<ndim){ // dim = nDim_field
                double position(0.);
                double min_limit(0.);
                double max_limit(0.);
                if (idim==0) {
                    position = cuParticles.position(0,iPart);
                    min_limit = min_local[0];
                    max_limit = max_local[0];
                }
                else { // look at square
                    position = cuParticles.position(1,iPart) * cuParticles.position(1,iPart)
                        + cuParticles.position(2,iPart) * cuParticles.position(2,iPart);
                    min_limit = min_local[0] * min_local[0];
                    max_limit = max_local[1] * max_local[1];
                }

                if ( position < min_limit ){
                    if ( neighbor_[idim][0]!=MPI_PROC_NULL) {
                        vecSpecies[ispec]->MPIbuff.part_index_send[idim][0].push_back( iPart );
                    }
                    //If particle is outside of the global domain (has no neighbor), it will not be put in a send buffer and will simply be deleted.
                    check = 1;
                }
                else if ( position >= max_limit ){
                    if( neighbor_[idim][1]!=MPI_PROC_NULL) {
                        vecSpecies[ispec]->MPIbuff.part_index_send[idim][1].push_back( iPart );
                    }
                    check = 1;
                }
                idim++;
            }
        }
    }

} // initExchParticles(... iDim)


// ---------------------------------------------------------------------------------------------------------------------
// For direction iDim, start exchange of number of particles
//   - vecPatch : used for intra-MPI process comm (direct copy using Particels::cp_particles)
//   - smpi     : inhereted from previous SmileiMPI::exchangeParticles()
// ---------------------------------------------------------------------------------------------------------------------
void Patch::initCommParticles(SmileiMPI* smpi, int ispec, Params& params, int iDim, VectorPatch * vecPatch)
{
    int h0 = (*vecPatch)(0)->hindex;
    /********************************************************************************/
    // Exchange number of particles to exchange to establish or not a communication
    /********************************************************************************/
    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
        if (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {
            vecSpecies[ispec]->MPIbuff.part_index_send_sz[iDim][iNeighbor] = (vecSpecies[ispec]->MPIbuff.part_index_send[iDim][iNeighbor]).size();

            if (is_a_MPI_neighbor(iDim, iNeighbor)) {
                //If neighbour is MPI ==> I send him the number of particles I'll send later.
                int tag = buildtag( hindex, iDim+1, iNeighbor+3 );
                      MPI_Isend( &(vecSpecies[ispec]->MPIbuff.part_index_send_sz[iDim][iNeighbor]), 1, MPI_INT, MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &(vecSpecies[ispec]->MPIbuff.srequest[iDim][iNeighbor]) );
            }
            else {
                //Else, I directly set the receive size to the correct value.
                (*vecPatch)( neighbor_[iDim][iNeighbor]- h0 )->vecSpecies[ispec]->MPIbuff.part_index_recv_sz[iDim][(iNeighbor+1)%2] = vecSpecies[ispec]->MPIbuff.part_index_send_sz[iDim][iNeighbor];
            }
        } // END of Send

        if (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
            if (is_a_MPI_neighbor(iDim, (iNeighbor+1)%2)) {
                //If other neighbour is MPI ==> I receive the number of particles I'll receive later.
                int tag = buildtag( neighbor_[iDim][(iNeighbor+1)%2], iDim+1, iNeighbor+3 );
                MPI_Irecv( &(vecSpecies[ispec]->MPIbuff.part_index_recv_sz[iDim][(iNeighbor+1)%2]), 1, MPI_INT, MPI_neighbor_[iDim][(iNeighbor+1)%2], tag, MPI_COMM_WORLD, &(vecSpecies[ispec]->MPIbuff.rrequest[iDim][(iNeighbor+1)%2]) );
            }
        }
    }//end loop on nb_neighbors.

} // initCommParticles(... iDim)


// ---------------------------------------------------------------------------------------------------------------------
// For direction iDim, finalize receive of number of particles and really send particles
//   - vecPatch : used for intra-MPI process comm (direct copy using Particels::cp_particles)
//   - smpi     : used smpi->periods_
// ---------------------------------------------------------------------------------------------------------------------
void Patch::CommParticles(SmileiMPI* smpi, int ispec, Params& params, int iDim, VectorPatch * vecPatch)
{
    Particles &cuParticles = (*vecSpecies[ispec]->particles);

    int n_part_send, n_part_recv;
    int h0 = (*vecPatch)(0)->hindex;
    double x_max = params.cell_length[iDim]*( params.n_space_global[iDim] );

    /********************************************************************************/
    // Wait for end of communications over number of particles
    /********************************************************************************/
    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
        MPI_Status sstat    [2];
        MPI_Status rstat    [2];
        if (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {
            if (is_a_MPI_neighbor(iDim, iNeighbor))
                MPI_Wait( &(vecSpecies[ispec]->MPIbuff.srequest[iDim][iNeighbor]), &(sstat[iNeighbor]) );
        }
        if (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
            if (is_a_MPI_neighbor(iDim, (iNeighbor+1)%2))  {
                MPI_Wait( &(vecSpecies[ispec]->MPIbuff.rrequest[iDim][(iNeighbor+1)%2]), &(rstat[(iNeighbor+1)%2]) );
                if (vecSpecies[ispec]->MPIbuff.part_index_recv_sz[iDim][(iNeighbor+1)%2]!=0) {
                    //If I receive particles over MPI, I initialize my receive buffer with the appropriate size.
                    vecSpecies[ispec]->MPIbuff.partRecv[iDim][(iNeighbor+1)%2].initialize( vecSpecies[ispec]->MPIbuff.part_index_recv_sz[iDim][(iNeighbor+1)%2], cuParticles);
                }
            }
        }
    }
    /********************************************************************************/
    // Proceed to effective Particles' communications
    /********************************************************************************/

    // Number of properties per particles = nDim_Particles + 3 + 1 + 1

    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {

        // n_part_send : number of particles to send to current neighbor
        n_part_send = (vecSpecies[ispec]->MPIbuff.part_index_send[iDim][iNeighbor]).size();
        if ( (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) && (n_part_send!=0) ) {
            // Enabled periodicity
            if (smpi->periods_[iDim]==1) {
                for (int iPart=0 ; iPart<n_part_send ; iPart++) {
                    if ( ( iNeighbor==0 ) &&  (Pcoordinates[iDim] == 0 ) &&( cuParticles.position(iDim,vecSpecies[ispec]->MPIbuff.part_index_send[iDim][iNeighbor][iPart]) < 0. ) ) {
                        cuParticles.position(iDim,vecSpecies[ispec]->MPIbuff.part_index_send[iDim][iNeighbor][iPart])     += x_max;
                    }
                    else if ( ( iNeighbor==1 ) &&  (Pcoordinates[iDim] == params.number_of_patches[iDim]-1 ) && ( cuParticles.position(iDim,vecSpecies[ispec]->MPIbuff.part_index_send[iDim][iNeighbor][iPart]) >= x_max ) ) {
                        cuParticles.position(iDim,vecSpecies[ispec]->MPIbuff.part_index_send[iDim][iNeighbor][iPart])     -= x_max;
                    }
                }
            }
            // Send particles
            if (is_a_MPI_neighbor(iDim, iNeighbor)) {
                // If MPI comm, first copy particles in the sendbuffer
                for (int iPart=0 ; iPart<n_part_send ; iPart++)
                    cuParticles.cp_particle(vecSpecies[ispec]->MPIbuff.part_index_send[iDim][iNeighbor][iPart], vecSpecies[ispec]->MPIbuff.partSend[iDim][iNeighbor]);
                // Then send particles
                int tag = buildtag( hindex, iDim+1, iNeighbor+3 );
                vecSpecies[ispec]->typePartSend[(iDim*2)+iNeighbor] = smpi->createMPIparticles( &(vecSpecies[ispec]->MPIbuff.partSend[iDim][iNeighbor]) );
                MPI_Isend( &((vecSpecies[ispec]->MPIbuff.partSend[iDim][iNeighbor]).position(0,0)), 1, vecSpecies[ispec]->typePartSend[(iDim*2)+iNeighbor], MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &(vecSpecies[ispec]->MPIbuff.srequest[iDim][iNeighbor]) );
            }
            else {
                //If not MPI comm, copy particles directly in the receive buffer
                for (int iPart=0 ; iPart<n_part_send ; iPart++)
                    cuParticles.cp_particle( vecSpecies[ispec]->MPIbuff.part_index_send[iDim][iNeighbor][iPart],((*vecPatch)( neighbor_[iDim][iNeighbor]- h0 )->vecSpecies[ispec]->MPIbuff.partRecv[iDim][(iNeighbor+1)%2]) );
            }
        } // END of Send

        n_part_recv = vecSpecies[ispec]->MPIbuff.part_index_recv_sz[iDim][(iNeighbor+1)%2];
        if ( (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
            if (is_a_MPI_neighbor(iDim, (iNeighbor+1)%2)) {
                // If MPI comm, receive particles in the recv buffer previously initialized.
                vecSpecies[ispec]->typePartRecv[(iDim*2)+iNeighbor] = smpi->createMPIparticles( &(vecSpecies[ispec]->MPIbuff.partRecv[iDim][(iNeighbor+1)%2]) );
                int tag = buildtag( neighbor_[iDim][(iNeighbor+1)%2], iDim+1 ,iNeighbor+3 );
                MPI_Irecv( &((vecSpecies[ispec]->MPIbuff.partRecv[iDim][(iNeighbor+1)%2]).position(0,0)), 1, vecSpecies[ispec]->typePartRecv[(iDim*2)+iNeighbor], MPI_neighbor_[iDim][(iNeighbor+1)%2], tag, MPI_COMM_WORLD, &(vecSpecies[ispec]->MPIbuff.rrequest[iDim][(iNeighbor+1)%2]) );
            }

        } // END of Recv

    } // END for iNeighbor

} // END CommParticles(... iDim)


// ---------------------------------------------------------------------------------------------------------------------
// For direction iDim, finalize receive of particles, temporary store particles if diagonalParticles
// And store recv particles at their definitive place.
// Call Patch::cleanup_sent_particles
//   - vecPatch : used for intra-MPI process comm (direct copy using Particels::cp_particles)
//   - smpi     : used smpi->periods_
// ---------------------------------------------------------------------------------------------------------------------
void Patch::finalizeCommParticles(SmileiMPI* smpi, int ispec, Params& params, int iDim, VectorPatch * vecPatch)
{

    int ndim = params.nDim_field;
    int idim, check;

    Particles &cuParticles = (*vecSpecies[ispec]->particles);

    int n_part_send, n_part_recv;

    /********************************************************************************/
    // Wait for end of communications over Particles
    /********************************************************************************/
    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
        MPI_Status sstat    [2];
        MPI_Status rstat    [2];

        n_part_send = vecSpecies[ispec]->MPIbuff.part_index_send[iDim][iNeighbor].size();
        n_part_recv = vecSpecies[ispec]->MPIbuff.part_index_recv_sz[iDim][(iNeighbor+1)%2];



        if ( (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) && (n_part_send!=0) ) {
            if (is_a_MPI_neighbor(iDim, iNeighbor)) {
                MPI_Wait( &(vecSpecies[ispec]->MPIbuff.srequest[iDim][iNeighbor]), &(sstat[iNeighbor]) );
                MPI_Type_free( &(vecSpecies[ispec]->typePartSend[(iDim*2)+iNeighbor]) );
            }
        }
        if ( (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
            if (is_a_MPI_neighbor(iDim, (iNeighbor+1)%2)) {
                MPI_Wait( &(vecSpecies[ispec]->MPIbuff.rrequest[iDim][(iNeighbor+1)%2]), &(rstat[(iNeighbor+1)%2]) );
                MPI_Type_free( &(vecSpecies[ispec]->typePartRecv[(iDim*2)+iNeighbor]) );
            }

            // Treat diagonalParticles
            if (iDim < ndim-1){ // No need to treat diag particles at last dimension.
                for (int iPart=n_part_recv-1 ; iPart>=0; iPart-- ) {
                    check = 0;
                    idim = iDim+1;//We check next dimension
                    while (check == 0 && idim<ndim){
                        //If particle not in the domain...
                        if ( (vecSpecies[ispec]->MPIbuff.partRecv[iDim][(iNeighbor+1)%2]).position(idim,iPart) < min_local[idim] ){
                            if (neighbor_[idim][0]!=MPI_PROC_NULL){ //if neighbour exists
                                //... copy it at the back of the local particle vector ...
                                (vecSpecies[ispec]->MPIbuff.partRecv[iDim][(iNeighbor+1)%2]).cp_particle(iPart, cuParticles);
                                //...adjust bmax or cell_keys ...
                                vecSpecies[ispec]->add_space_for_a_particle();
                                //... and add its index to the particles to be sent later...
                                vecSpecies[ispec]->MPIbuff.part_index_send[idim][0].push_back( cuParticles.size()-1 );
                                //..without forgeting to add it to the list of particles to clean.
                                vecSpecies[ispec]->addPartInExchList(cuParticles.size()-1);
                            }
                            //Remove it from receive buffer.
                            (vecSpecies[ispec]->MPIbuff.partRecv[iDim][(iNeighbor+1)%2]).erase_particle(iPart);
                            vecSpecies[ispec]->MPIbuff.part_index_recv_sz[iDim][(iNeighbor+1)%2]--;
                            check = 1;
                        }
                        //Other side of idim
                        else if ( (vecSpecies[ispec]->MPIbuff.partRecv[iDim][(iNeighbor+1)%2]).position(idim,iPart) >= max_local[idim]) {
                            if (neighbor_[idim][1]!=MPI_PROC_NULL){ //if neighbour exists
                                (vecSpecies[ispec]->MPIbuff.partRecv[iDim][(iNeighbor+1)%2]).cp_particle(iPart, cuParticles);
                                //...adjust bmax or cell_keys ...
                                vecSpecies[ispec]->add_space_for_a_particle();
                                vecSpecies[ispec]->MPIbuff.part_index_send[idim][1].push_back( cuParticles.size()-1 );
                                vecSpecies[ispec]->addPartInExchList(cuParticles.size()-1);
                            }
                            (vecSpecies[ispec]->MPIbuff.partRecv[iDim][(iNeighbor+1)%2]).erase_particle(iPart);
                            vecSpecies[ispec]->MPIbuff.part_index_recv_sz[iDim][(iNeighbor+1)%2]--;
                            check = 1;
                        }
                        idim++;
                    }
                }
            }//If not last dim for diagonal particles.

        } //If received something
    } //loop i Neighbor

    //La recopie finale doit se faire au traitement de la dernière dimension seulement !!
    if (iDim == ndim-1){

        //vecSpecies[ispec]->sort_part(params);

        // For DynamicV
        if (vecSpecies[ispec]->vectorized_operators)
        {
            vecSpecies[ispec]->sort_part(params);
        }
        else
        {
            vecSpecies[ispec]->Species::sort_part(params);
        }

    }//End Recv_buffers ==> particles


} // finalizeCommParticles(... iDim)


void Patch::cleanParticlesOverhead(Params& params)
{
    int ndim = params.nDim_field;
    for (unsigned int ispec=0 ; ispec<vecSpecies.size() ; ispec++) {
        Particles &cuParticles = (*vecSpecies[ispec]->particles);

        for (int idim = 0; idim < ndim; idim++){
            for ( int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++ ) {
                vecSpecies[ispec]->MPIbuff.partRecv[idim][iNeighbor].clear();
                vecSpecies[ispec]->MPIbuff.partRecv[idim][iNeighbor].shrink_to_fit(ndim);
                vecSpecies[ispec]->MPIbuff.partSend[idim][iNeighbor].clear();
                vecSpecies[ispec]->MPIbuff.partSend[idim][iNeighbor].shrink_to_fit(ndim);
                vecSpecies[ispec]->MPIbuff.part_index_send[idim][iNeighbor].clear();
                vector<int>(vecSpecies[ispec]->MPIbuff.part_index_send[idim][iNeighbor]).swap(vecSpecies[ispec]->MPIbuff.part_index_send[idim][iNeighbor]);
            }
        }

        cuParticles.shrink_to_fit(ndim);
    }

}

// ---------------------------------------------------------------------------------------------------------------------
// Clear vecSpecies[]->indexes_of_particles_to_exchange, suppress particles send and manage memory
// ---------------------------------------------------------------------------------------------------------------------
void Patch::cleanup_sent_particles(int ispec, std::vector<int>* indexes_of_particles_to_exchange)
{
    /********************************************************************************/
    // Delete Particles included in the index of particles to exchange. Assumes indexes are sorted.
    /********************************************************************************/
    int ii, iPart;
    std::vector<int>* cubmin = &vecSpecies[ispec]->bmin;
    std::vector<int>* cubmax = &vecSpecies[ispec]->bmax;
    Particles &cuParticles = (*vecSpecies[ispec]->particles);


    // Push lost particles at the end of bins
    for (unsigned int ibin = 0 ; ibin < (*cubmax).size() ; ibin++ ) {
        ii = (*indexes_of_particles_to_exchange).size()-1;
        if (ii >= 0) { // Push lost particles to the end of the bin
            iPart = (*indexes_of_particles_to_exchange)[ii];
            while (iPart >= (*cubmax)[ibin] && ii > 0) {
                ii--;
                iPart = (*indexes_of_particles_to_exchange)[ii];
            }
            while (iPart == (*cubmax)[ibin]-1 && iPart >= (*cubmin)[ibin] && ii > 0) {
                (*cubmax)[ibin]--;
                ii--;
                iPart = (*indexes_of_particles_to_exchange)[ii];
            }
            while (iPart >= (*cubmin)[ibin] && ii > 0) {
                cuParticles.overwrite_part((*cubmax)[ibin]-1, iPart );
                (*cubmax)[ibin]--;
                ii--;
                iPart = (*indexes_of_particles_to_exchange)[ii];
            }
            if (iPart >= (*cubmin)[ibin] && iPart < (*cubmax)[ibin]) { //On traite la dernière particule (qui peut aussi etre la premiere)
                cuParticles.overwrite_part((*cubmax)[ibin]-1, iPart );
                (*cubmax)[ibin]--;
            }
        }
    }


    //Shift the bins in memory
    //Warning: this loop must be executed sequentially. Do not use openMP here.
    for (int unsigned ibin = 1 ; ibin < (*cubmax).size() ; ibin++ ) { //First bin don't need to be shifted
        ii = (*cubmin)[ibin]-(*cubmax)[ibin-1]; // Shift the bin in memory by ii slots.
        iPart = min(ii,(*cubmax)[ibin]-(*cubmin)[ibin]); // Number of particles we have to shift = min (Nshift, Nparticle in the bin)
        if(iPart > 0) cuParticles.overwrite_part((*cubmax)[ibin]-iPart,(*cubmax)[ibin-1],iPart);
        (*cubmax)[ibin] -= ii;
        (*cubmin)[ibin] = (*cubmax)[ibin-1];
    }

} // END cleanup_sent_particles

void Patch::testSumField( Field* field, int iDim )
{
    MPI_Status sstat    [2];
    MPI_Status rstat    [2];
    int tagp(0);
    if (field->name == "Jx") tagp = 1;
    if (field->name == "Jy") tagp = 2;
    if (field->name == "Jz") tagp = 3;
    if (field->name == "Rho") tagp = 4;

    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
        if ( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
            int received(false);
            MPI_Test( &(field->MPIbuff.srequest[iDim][iNeighbor]), &received, &(sstat[iNeighbor]) );
        }
        if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {
            int received(false);

            int tag  = buildtag( neighbor_[iDim][(iNeighbor+1)%2], iDim, iNeighbor, tagp );
            MPI_Iprobe(MPI_neighbor_[iDim][(iNeighbor+1)%2], tag, MPI_COMM_WORLD, &received, &(rstat[(iNeighbor+1)%2]));
            //MPI_Test( &(field->MPIbuff.rrequest[iDim][(iNeighbor+1)%2]), &received, &(rstat[(iNeighbor+1)%2]) );
        }
    }

} // END testSumField
