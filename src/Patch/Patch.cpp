
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
#include "InterpolatorFactory.h"
#include "ProjectorFactory.h"
#include "DiagnosticFactory.h"

using namespace std;


// ---------------------------------------------------------------------------------------------------------------------
// Patch constructor :
//   Called by PatchXD constructor which will finalize initialization
// ---------------------------------------------------------------------------------------------------------------------
Patch::Patch(Params& params, SmileiMPI* smpi, unsigned int ipatch, unsigned int n_moved)
{
    
    hindex = ipatch;
    nDim_fields_ = params.nDim_field;
    
    initStep1(params);
    
} // END Patch::Patch




// Cloning patch constructor
Patch::Patch(Patch* patch, Params& params, SmileiMPI* smpi, unsigned int ipatch, unsigned int n_moved, bool with_particles = true) {
    
    hindex = ipatch;
    nDim_fields_ = patch->nDim_fields_;
    
    initStep1(params);
    
}


void Patch::initStep1(Params& params)
{
    // for nDim_fields = 1 : bug if Pcoordinates.size = 1 !! 
    //Pcoordinates.resize(nDim_fields_);
    Pcoordinates.resize( 2 );
    
    // else if ( params.geometry == "2d3v" ) {
    //     Pcoordinates.resize(3);
    //     generalhilbertindexinv(params.mi[0], params.mi[1], params.mi[2], &Pcoordinates[0], &Pcoordinates[1], &Pcoordinates[2], hindex);
    // }
    
    nbNeighbors_ = 2;
    neighbor_.resize(nDim_fields_);
    corner_neighbor_.resize(params.nDim_field);
    for ( int iDim = 0 ; iDim < nDim_fields_ ; iDim++ ) {
        neighbor_[iDim].resize(2,MPI_PROC_NULL);
        corner_neighbor_[iDim].resize(2,MPI_PROC_NULL);
    }
    MPI_neighbor_.resize(nDim_fields_);
    for ( int iDim = 0 ; iDim < nDim_fields_; iDim++ ) {
        MPI_neighbor_[iDim].resize(2,MPI_PROC_NULL);
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
    cell_starting_global_index.resize(params.nDim_field, 0);
    for (unsigned int i = 0 ; i<params.nDim_field ; i++) {
        min_local[i] =  Pcoordinates[i]   *params.n_space[i]*params.cell_length[i];
        max_local[i] = (Pcoordinates[i]+1)*params.n_space[i]*params.cell_length[i];
        cell_starting_global_index[i] += Pcoordinates[i]*params.n_space[i];
        cell_starting_global_index[i] -= params.oversize[i];
    }
    
    cell_starting_global_index[0] += n_moved;
    min_local[0] += n_moved*params.cell_length[0];
    max_local[0] += n_moved*params.cell_length[0];
}


void Patch::finishCreation( Params& params, SmileiMPI* smpi ) {
    // initialize vector of Species (virtual)
    vecSpecies = SpeciesFactory::createVector(params, this);
    
    // initialize the electromagnetic fields (virtual)
    EMfields   = ElectroMagnFactory::create(params, vecSpecies, this);
    
    // interpolation operator (virtual)
    Interp     = InterpolatorFactory::create(params, this); // + patchId -> idx_domain_begin (now = ref smpi)
    // projection operator (virtual)
    Proj       = ProjectorFactory::create(params, this);    // + patchId -> idx_domain_begin (now = ref smpi)
    
    // Initialize the collisions
    vecCollisions = Collisions::create(params, this, vecSpecies);
    
    // Initialize the particle walls
    partWalls = new PartWalls(params, this);
    
    // Initialize the probes
    probes = DiagnosticFactory::createProbes();
    
    createType(params);
}


void Patch::finishCloning( Patch* patch, Params& params, SmileiMPI* smpi, bool with_particles = true ) {
    // clone vector of Species (virtual)
    vecSpecies = SpeciesFactory::cloneVector(patch->vecSpecies, params, this, with_particles);
    
    // clone the electromagnetic fields (virtual)
    EMfields   = ElectroMagnFactory::clone(patch->EMfields, params, vecSpecies, this);
    
    // interpolation operator (virtual)
    Interp     = InterpolatorFactory::create(params, this);
    // projection operator (virtual)
    Proj       = ProjectorFactory::create(params, this);
    
    // clone the collisions
    vecCollisions = Collisions::clone(patch->vecCollisions, params);
    
    // clone the particle walls
    partWalls = new PartWalls(patch->partWalls, this);
    
    // clone the probes
    probes = DiagnosticFactory::cloneProbes(patch->probes);
    
    createType(params);
}


// ---------------------------------------------------------------------------------------------------------------------
// Delete Patch members
// ---------------------------------------------------------------------------------------------------------------------
Patch::~Patch() {
    
    for(unsigned int i=0; i<vecCollisions.size(); i++) delete vecCollisions[i];
    vecCollisions.clear();
    
    delete partWalls;
    
    delete Proj;
    delete Interp;
    
    delete EMfields;
    for (unsigned int ispec=0 ; ispec<vecSpecies.size(); ispec++) delete vecSpecies[ispec];
    vecSpecies.clear();
    
} // END Patch::~Patch


// ---------------------------------------------------------------------------------------------------------------------
// Compute MPI rank of patch neigbors and current patch
// ---------------------------------------------------------------------------------------------------------------------
void Patch::updateMPIenv(SmileiMPI* smpi)
{
    MPI_me_ = smpi->smilei_rk;
    
    for (int iDim = 0 ; iDim < nDim_fields_ ; iDim++)
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++)
            MPI_neighbor_[iDim][iNeighbor] = smpi->hrank(neighbor_[iDim][iNeighbor]);
    
//        cout << "\n\tPatch Corner decomp : " << corner_neighbor_[0][1] << "\t" << neighbor_[1][1]  << "\t" << corner_neighbor_[1][1] << endl;
//        cout << "\tPatch Corner decomp : " << neighbor_[0][0] << "\t" << hindex << "\t" << neighbor_[0][1] << endl;
//        cout << "\tPatch Corner decomp : " << corner_neighbor_[0][0] << "\t" << neighbor_[1][0]  << "\t" << corner_neighbor_[1][0] << endl;
//        

//        cout << "\n\tMPI Corner decomp : " << "MPI_PROC_NULL" << "\t" << MPI_neighbor_[2][1]  << "\t" << "MPI_PROC_NULL" << endl << endl;
//
//        cout << "\n\tMPI Corner decomp : " << "MPI_PROC_NULL" << "\t" << MPI_neighbor_[1][1]  << "\t" << "MPI_PROC_NULL" << endl;
//        cout << "\tMPI Corner decomp : " << MPI_neighbor_[0][0] << "\t" << smpi->getRank() << "\t" << MPI_neighbor_[0][1] << endl;
//        cout << "\tMPI Corner decomp : " << "MPI_PROC_NULL" << "\t" << MPI_neighbor_[1][0]  << "\t" << "MPI_PROC_NULL" << endl << endl;
//
//        cout << "\n\tMPI Corner decomp : " << "MPI_PROC_NULL" << "\t" << MPI_neighbor_[2][0]  << "\t" << "MPI_PROC_NULL" << endl;
    
} // END updateMPIenv


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
    double xmax[3]; 
    
    for (int iDim=0 ; iDim < ndim ; iDim++){
        xmax[iDim] = params.cell_length[iDim]*( params.n_space_global[iDim] );
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
            vecSpecies[ispec]->MPIbuff.partRecv[iDim][iNeighbor].initialize(0,cuParticles);
            vecSpecies[ispec]->MPIbuff.partSend[iDim][iNeighbor].initialize(0,cuParticles);
            vecSpecies[ispec]->MPIbuff.part_index_send[iDim][iNeighbor].resize(0);
            vecSpecies[ispec]->MPIbuff.part_index_recv_sz[iDim][iNeighbor] = 0;
        }
    }
 
    int n_part_send = (*indexes_of_particles_to_exchange).size();
    //int n_part_recv;
        
    int iPart;

    // Define where particles are going 
    //Put particles in the send buffer it belongs to. Priority to lower dimensions.
    for (int i=0 ; i<n_part_send ; i++) {
        iPart = (*indexes_of_particles_to_exchange)[i];
        check = 0;
        idim = 0;
        //Put indexes of particles in the first direction they will be exchanged and correct their position according to periodicity for the first exchange only.
        while (check == 0 && idim<ndim){
            if ( cuParticles.position(idim,iPart) < min_local[idim]){
                if ( neighbor_[idim][0]!=MPI_PROC_NULL) { 
                    vecSpecies[ispec]->MPIbuff.part_index_send[idim][0].push_back( iPart );
                    if (smpi->periods_[idim]==1 && Pcoordinates[idim] == 0) {
                        cuParticles.position(idim,iPart)     += xmax[idim];
                    }
                }
                //If particle is outside of the global domain (has no neighbor), it will not be put in a send buffer and will simply be deleted.
                check = 1;
            }
            else if ( cuParticles.position(idim,iPart) >= max_local[idim]){
                if( neighbor_[idim][1]!=MPI_PROC_NULL) { 
                    vecSpecies[ispec]->MPIbuff.part_index_send[idim][1].push_back( iPart );
                    if (smpi->periods_[idim]==1 && (int)Pcoordinates[idim] == params.number_of_patches[idim]-1) {
                        cuParticles.position(idim,iPart)     -= xmax[idim];
                    }
                }
                check = 1;
            }
            idim++;
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
    MPI_Datatype typePartSend, typePartRecv;
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
                    else if ( ( iNeighbor==1 ) &&  ((int)Pcoordinates[iDim] == params.number_of_patches[iDim]-1 ) && ( cuParticles.position(iDim,vecSpecies[ispec]->MPIbuff.part_index_send[iDim][iNeighbor][iPart]) >= x_max ) ) {
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
                typePartSend = smpi->createMPIparticles( &(vecSpecies[ispec]->MPIbuff.partSend[iDim][iNeighbor]) );
                MPI_Isend( &((vecSpecies[ispec]->MPIbuff.partSend[iDim][iNeighbor]).position(0,0)), 1, typePartSend, MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &(vecSpecies[ispec]->MPIbuff.srequest[iDim][iNeighbor]) );
                MPI_Type_free( &typePartSend );
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
                typePartRecv = smpi->createMPIparticles( &(vecSpecies[ispec]->MPIbuff.partRecv[iDim][(iNeighbor+1)%2]) );
                int tag = buildtag( neighbor_[iDim][(iNeighbor+1)%2], iDim+1 ,iNeighbor+3 );
                MPI_Irecv( &((vecSpecies[ispec]->MPIbuff.partRecv[iDim][(iNeighbor+1)%2]).position(0,0)), 1, typePartRecv, MPI_neighbor_[iDim][(iNeighbor+1)%2], tag, MPI_COMM_WORLD, &(vecSpecies[ispec]->MPIbuff.rrequest[iDim][(iNeighbor+1)%2]) );
                MPI_Type_free( &typePartRecv );
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
    double xmax[3]; 
    for (int idim=0 ; idim < ndim ; idim++)
        xmax[idim] = params.cell_length[idim]*( params.n_space_global[idim] );


    std::vector<int>* indexes_of_particles_to_exchange = &vecSpecies[ispec]->indexes_of_particles_to_exchange;

    std::vector<int>* cubmin = &vecSpecies[ispec]->bmin;
    std::vector<int>* cubmax = &vecSpecies[ispec]->bmax;

    int nmove,lmove,ii; // local, OK
    int shift[(*cubmax).size()+1];//how much we need to shift each bin in order to leave room for the new particles
    double dbin;
        
    dbin = params.cell_length[0]*params.clrw; //width of a bin.
    for (unsigned int j=0; j<(*cubmax).size()+1 ;j++){
      shift[j]=0;
    }

    int n_part_send, n_part_recv, n_particles;

    /********************************************************************************/
    // Wait for end of communications over Particles
    /********************************************************************************/
    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
        MPI_Status sstat    [2];
        MPI_Status rstat    [2];
                
        n_part_send = vecSpecies[ispec]->MPIbuff.part_index_send[iDim][iNeighbor].size();
        n_part_recv = vecSpecies[ispec]->MPIbuff.part_index_recv_sz[iDim][(iNeighbor+1)%2];
               

 
        if ( (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) && (n_part_send!=0) ) {
            if (is_a_MPI_neighbor(iDim, iNeighbor))
                MPI_Wait( &(vecSpecies[ispec]->MPIbuff.srequest[iDim][iNeighbor]), &(sstat[iNeighbor]) );
        }
        if ( (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
            if (is_a_MPI_neighbor(iDim, (iNeighbor+1)%2))
                MPI_Wait( &(vecSpecies[ispec]->MPIbuff.rrequest[iDim][(iNeighbor+1)%2]), &(rstat[(iNeighbor+1)%2]) );     

            // Treat diagonalParticles
            if (iDim < ndim-1){ // No need to treat diag particles at last dimension.
                for (int iPart=n_part_recv-1 ; iPart>=0; iPart-- ) {
                    check = 0;
                    idim = iDim+1;//We check next dimension
                    while (check == 0 && idim<ndim){
                        //If particle not in the domain...
                        if ( (vecSpecies[ispec]->MPIbuff.partRecv[iDim][(iNeighbor+1)%2]).position(idim,iPart) < min_local[idim] ){  
                            if (neighbor_[idim][0]!=MPI_PROC_NULL){ //if neighbour exists
                                //...Deal with periodicity...
                                    if (smpi->periods_[idim]==1 && Pcoordinates[idim] == 0) {
                                    (vecSpecies[ispec]->MPIbuff.partRecv[iDim][(iNeighbor+1)%2]).position(idim,iPart)     += xmax[idim];
                                }
                                //... copy it at the back of the local particle vector ...
                                (vecSpecies[ispec]->MPIbuff.partRecv[iDim][(iNeighbor+1)%2]).cp_particle(iPart, cuParticles);
                                //...adjust bmax ...
                                (*cubmax)[(*cubmax).size()-1]++;
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
                                    if (smpi->periods_[idim]==1 && (int)Pcoordinates[idim] == params.number_of_patches[idim]-1) {
                                    (vecSpecies[ispec]->MPIbuff.partRecv[iDim][(iNeighbor+1)%2]).position(idim,iPart)     -= xmax[idim];
                                }
                                (vecSpecies[ispec]->MPIbuff.partRecv[iDim][(iNeighbor+1)%2]).cp_particle(iPart, cuParticles);
                                (*cubmax)[(*cubmax).size()-1]++;
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
   
        //We have stored in indexes_of_particles_to_exchange the list of all particles that needs to be removed.
        cleanup_sent_particles(ispec, indexes_of_particles_to_exchange);
        (*indexes_of_particles_to_exchange).clear();
        cuParticles.erase_particle_trail((*cubmax).back());

        //Evaluation of the necessary shift of all bins.
        for (unsigned int j=0; j<(*cubmax).size()+1 ;j++){
            shift[j]=0;
        }

        //idim=0
        shift[1] += vecSpecies[ispec]->MPIbuff.part_index_recv_sz[0][0];//Particles coming from south all go to bin 0 and shift all the other bins.
        shift[(*cubmax).size()] += vecSpecies[ispec]->MPIbuff.part_index_recv_sz[0][1];//Used only to count the total number of particles arrived.
        //idim>0
        for (idim = 1; idim < ndim; idim++){
            for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
                n_part_recv = vecSpecies[ispec]->MPIbuff.part_index_recv_sz[idim][iNeighbor];
                for (unsigned int j=0; j<(unsigned int)n_part_recv ;j++){
                    //We first evaluate how many particles arrive in each bin.
                        ii = int((vecSpecies[ispec]->MPIbuff.partRecv[idim][iNeighbor].position(0,j)-min_local[0])/dbin);//bin in which the particle goes.
                    shift[ii+1]++; // It makes the next bins shift.
                }
            }
        }


        //Must be done sequentially
        for (unsigned int j=1; j<(*cubmax).size()+1;j++){ //bin 0 is not shifted.Last element of shift stores total number of arriving particles.
            shift[j]+=shift[j-1];
        }
        //Make room for new particles
        if (shift[(*cubmax).size()]) {
          //! vecor::resize of Charge crashed ! Temporay solution : push_back / Particle
          //cuParticles.initialize( cuParticles.size()+shift[(*cubmax).size()], cuParticles.Position.size() );
          for (int inewpart=0 ; inewpart<shift[(*cubmax).size()] ; inewpart++) cuParticles.create_particle();
        }
            
        //Shift bins, must be done sequentially
        for (unsigned int j=(*cubmax).size()-1; j>=1; j--){
            n_particles = (*cubmax)[j]-(*cubmin)[j]; //Nbr of particle in this bin
            nmove = min(n_particles,shift[j]); //Nbr of particles to move
            lmove = max(n_particles,shift[j]); //How far particles must be shifted
            if (nmove>0) cuParticles.overwrite_part((*cubmin)[j], (*cubmin)[j]+lmove, nmove);
            (*cubmin)[j] += shift[j];
            (*cubmax)[j] += shift[j];
        }
            
        //Space has been made now to write the arriving particles into the correct bins
        //idim == 0  is the easy case, when particles arrive either in first or last bin.
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
            n_part_recv = vecSpecies[ispec]->MPIbuff.part_index_recv_sz[0][iNeighbor];
            if ( (neighbor_[0][iNeighbor]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
                ii = iNeighbor*((*cubmax).size()-1);//0 if iNeighbor=0(particles coming from West) and (*cubmax).size()-1 otherwise.
                vecSpecies[ispec]->MPIbuff.partRecv[0][iNeighbor].overwrite_part(0, cuParticles,(*cubmax)[ii],n_part_recv);
                (*cubmax)[ii] += n_part_recv ;
            }
        }
        //idim > 0; this is the difficult case, when particles can arrive in any bin.
        for (idim = 1; idim < ndim; idim++){
            for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
                n_part_recv = vecSpecies[ispec]->MPIbuff.part_index_recv_sz[idim][iNeighbor];
                if ( (neighbor_[idim][iNeighbor]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
                        for(unsigned int j=0; j<(unsigned int)n_part_recv; j++){
                            ii = int((vecSpecies[ispec]->MPIbuff.partRecv[iDim][iNeighbor].position(0,j)-min_local[0])/dbin);//bin in which the particle goes.
                            vecSpecies[ispec]->MPIbuff.partRecv[idim][iNeighbor].overwrite_part(j, cuParticles,(*cubmax)[ii]);
                            (*cubmax)[ii] ++ ;
                        }
                }
            }
        }
    }//End Recv_buffers ==> particles

} // finalizeCommParticles(... iDim)


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
