
#include "Patch.h"

#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "Hilbert_functions.h"
#include "PatchesFactory.h"
#include "Species.h"
#include "Particles.h"
#include "SmileiIOFactory.h"

using namespace std;

//int buildtag(int send, int recv);

Patch::Patch(PicParams& params, DiagParams &diag_params, LaserParams& laser_params, SmileiMPI* smpi, unsigned int ipatch, unsigned int n_moved) {


//Neighborhood definition in 2D:
//   
//   Y axis
//   ^
//   |   6       7          8
//   |   3    4(self)       5
//       0       1          2    --> X axis

        int xcall, ycall;
        hindex = ipatch;
#ifdef _DEBUGPATCH
	std::cout << smpi->getRank() << ", mypatch is : " << hindex << std::endl;        
#endif
        if ( params.geometry == "1d3v" ) {
            Pcoordinates.resize(1);
            Pcoordinates[0] = hindex;
	    MPI_neighborhood_.resize(3);
	    patch_neighborhood_.resize(3);
        }
        else if ( params.geometry == "2d3v" ) {
            Pcoordinates.resize(2);
            generalhilbertindexinv(params.mi[0], params.mi[1], &Pcoordinates[0], &Pcoordinates[1], hindex);

	    /////////////////////////////////////
	    //  Define local domain
	    /////////////////////////////////////

	    MPI_neighborhood_.resize(9);
	    patch_neighborhood_.resize(9);
        }
        else {
            Pcoordinates.resize(3);
            generalhilbertindexinv(params.mi[0], params.mi[1], params.mi[2], &Pcoordinates[0], &Pcoordinates[1], &Pcoordinates[2], hindex);
	    MPI_neighborhood_.resize(27);
	    patch_neighborhood_.resize(27);
        }

	//std::cout << "Coordonnées de " << ipatch << " : " << Pcoordinates[0] << " " << Pcoordinates[1] << std::endl;
	nbNeighbors_ = 2;
	neighbor_.resize(params.nDim_field);
	corner_neighbor_.resize(params.nDim_field);
	for ( int iDim = 0 ; iDim < params.nDim_field ; iDim++ ) {
	    neighbor_[iDim].resize(2,MPI_PROC_NULL);
	    corner_neighbor_[iDim].resize(2,MPI_PROC_NULL);
	}
	MPI_neighbor_.resize(params.nDim_field);
	MPI_corner_neighbor_.resize(params.nDim_field);
	for ( int iDim = 0 ; iDim < params.nDim_field ; iDim++ ) {
	    MPI_neighbor_[iDim].resize(2,MPI_PROC_NULL);
	    MPI_corner_neighbor_[iDim].resize(2,MPI_PROC_NULL);
	}

        xcall = Pcoordinates[0]-1;
        ycall = Pcoordinates[1];
	if (params.bc_em_type_long=="periodic" && xcall < 0) xcall += (1<<params.mi[0]);
	neighbor_[0][0] = generalhilbertindex( params.mi[0], params.mi[1], xcall, ycall);
#ifdef _PATCH_DEBUG
        cout << xcall << " " << ycall << " " << neighbor_[0][0] << endl;
#endif
        xcall = Pcoordinates[0]+1;
	if (params.bc_em_type_long=="periodic" && xcall >= (1<<params.mi[0])) xcall -= (1<<params.mi[0]);
	neighbor_[0][1] = generalhilbertindex( params.mi[0], params.mi[1], xcall, ycall);
        xcall = Pcoordinates[0];
        ycall = Pcoordinates[1]-1;
	if (params.bc_em_type_trans=="periodic" && ycall < 0) ycall += (1<<params.mi[1]);
	neighbor_[1][0] = generalhilbertindex( params.mi[0], params.mi[1], xcall, ycall);
        ycall = Pcoordinates[1]+1;
	if (params.bc_em_type_trans=="periodic" && ycall >= (1<<params.mi[1])) ycall -= (1<<params.mi[1]);
	neighbor_[1][1] = generalhilbertindex( params.mi[0], params.mi[1], xcall, ycall);


        xcall = Pcoordinates[0]+1;
	if (params.bc_em_type_long=="periodic" && xcall >= (1<<params.mi[0])) xcall -= (1<<params.mi[0]);
	corner_neighbor_[1][1] = generalhilbertindex( params.mi[0], params.mi[1], xcall, ycall);
        xcall = Pcoordinates[0]-1;
	if (params.bc_em_type_long=="periodic" && xcall < 0) xcall += (1<<params.mi[0]);
	corner_neighbor_[0][1] = generalhilbertindex( params.mi[0], params.mi[1], xcall, ycall);
        ycall = Pcoordinates[1]-1;
	if (params.bc_em_type_trans=="periodic" && ycall < 0) ycall += (1<<params.mi[1]);
	corner_neighbor_[0][0] = generalhilbertindex( params.mi[0], params.mi[1], xcall, ycall);
        xcall = Pcoordinates[0]+1;
	if (params.bc_em_type_long=="periodic" && xcall >= (1<<params.mi[0])) xcall -= (1<<params.mi[0]);
	corner_neighbor_[1][0] = generalhilbertindex( params.mi[0], params.mi[1], xcall, ycall);


        patch_neighborhood_[0] = corner_neighbor_[0][0];
        patch_neighborhood_[1] = neighbor_[1][0];
        patch_neighborhood_[2] = corner_neighbor_[1][0];
        patch_neighborhood_[3] = neighbor_[0][0];
        patch_neighborhood_[4] = hindex;
        patch_neighborhood_[5] = neighbor_[0][1];
        patch_neighborhood_[6] = corner_neighbor_[0][1];
        patch_neighborhood_[7] = neighbor_[1][1];
        patch_neighborhood_[8] = corner_neighbor_[1][1];

	updateMPIenv(smpi);

	createType(params);

	
	//std::cout << "Voisin dir 0 : " << ipatch << " : " <<  neighbor_[0][0] << " " <<  neighbor_[0][1] << std::endl;
	//std::cout << "Voisin dir 1 : " << ipatch << " : " <<  neighbor_[1][0] << " " <<  neighbor_[1][1] << std::endl;

	min_local.resize(params.nDim_field, 0.);
	max_local.resize(params.nDim_field, 0.);
	cell_starting_global_index.resize(params.nDim_field, 0);
	for (int i = 0 ; i<params.nDim_field ; i++) {
	    min_local[i] = Pcoordinates[i]*params.n_space[i]*params.cell_length[i];
	    max_local[i] = min_local[i] + params.n_space[i]*params.cell_length[i];
	    cell_starting_global_index[i] += Pcoordinates[i]*params.n_space[i];
	    cell_starting_global_index[i] -= params.oversize[i];
	}

	cell_starting_global_index[0] += n_moved;
	min_local[0] += n_moved*params.cell_length[0];
	max_local[0] += n_moved*params.cell_length[0];

	vecSpecies = SpeciesFactory::createVector(params, this);

	/* // + min_loc/cell_index(ref smpi,  & sort) // OK through this */
	// + new n_space -> in PatchFactory
	// patchID : ok through coord
	// create Pos : OK

	// -> partBoundCond : min/max_loc (smpi)
	EMfields   = ElectroMagnFactory::create(params, laser_params, this);
	// + patchId + new n_space (now = params by smpi) + BC
	// -> Neighbors to define !!
	
	Interp     = InterpolatorFactory::create(params, this);               // + patchId -> idx_domain_begin (now = ref smpi)
	Proj       = ProjectorFactory::create(params, this);                  // + patchId -> idx_domain_begin (now = ref smpi)

	Diags = new Diagnostic(params,diag_params, this);

	sio = SmileiIOFactory::create(params, diag_params, this);
	
};

void Patch::updateMPIenv(SmileiMPI* smpi)
{
	for ( int z = 0 ; z < 1+2*(2 == 3) ; z++ ) {
	    for ( int y = 0 ; y < 1+2*(2 >= 2) ; y++ ) {
	        for ( int x = 0 ; x < 3 ; x++ ) {
		    MPI_neighborhood_[z*9+y*3+x] = smpi->hrank(patch_neighborhood_[z*9+y*3+x]);
                }
            }
        }

#ifdef _PATCH_DEBUG
	cout << "\n\tPatch Corner decomp : " << corner_neighbor_[0][1] << "\t" << neighbor_[1][1]  << "\t" << corner_neighbor_[1][1] << endl;
	cout << "\tPatch Corner decomp : " << neighbor_[0][0] << "\t" << hindex << "\t" << neighbor_[0][1] << endl;
	cout << "\tPatch Corner decomp : " << corner_neighbor_[0][0] << "\t" << neighbor_[1][0]  << "\t" << corner_neighbor_[1][0] << endl;

	cout << "\n\tMPI Corner decomp : " << MPI_neighborhood_[6] << "\t" << MPI_neighborhood_[7]  << "\t" << MPI_neighborhood_[8] << endl;
	cout << "\tMPI Corner decomp : " << MPI_neighborhood_[3] << "\t" << MPI_neighborhood_[4] << "\t" << MPI_neighborhood_[5] << endl;
	cout << "\tMPI Corner decomp : " << MPI_neighborhood_[0] << "\t" << MPI_neighborhood_[1]  << "\t" << MPI_neighborhood_[2] << endl;
#endif

	// Redundant temporary solution, to introduce, MPI-Patched features
	MPI_neighbor_[0][0] = MPI_neighborhood_[3];
	MPI_neighbor_[0][1] = MPI_neighborhood_[5];
	MPI_neighbor_[1][0] = MPI_neighborhood_[1];
	MPI_neighbor_[1][1] = MPI_neighborhood_[7];
	MPI_corner_neighbor_[0][0] = MPI_neighborhood_[0];
	MPI_corner_neighbor_[0][1] = MPI_neighborhood_[6];
	MPI_corner_neighbor_[1][0] = MPI_neighborhood_[2];
	MPI_corner_neighbor_[1][1] = MPI_neighborhood_[8];

#ifdef _PATCH_DEBUG
	cout << "\n\tMPI Corner decomp : " << MPI_corner_neighbor_[0][1] << "\t" << MPI_neighbor_[1][1]  << "\t" << MPI_corner_neighbor_[1][1] << endl;
	cout << "\tMPI Corner decomp : " << MPI_neighbor_[0][0] << "\t" << smpi->getRank() << "\t" << MPI_neighbor_[0][1] << endl;
	cout << "\tMPI Corner decomp : " << MPI_corner_neighbor_[0][0] << "\t" << MPI_neighbor_[1][0]  << "\t" << MPI_corner_neighbor_[1][0] << endl;
#endif

}


void Patch::dynamics(double time_dual, PicParams &params, SimWindow* simWindow, int diag_flag)
{
    for (unsigned int ispec=0 ; ispec<params.n_species; ispec++) {
	if ( vecSpecies[ispec]->isProj(time_dual, simWindow) || diag_flag  ){    
	    vecSpecies[ispec]->dynamics(time_dual, ispec, EMfields, Interp, Proj, params, diag_flag);
	}
    }

}

//void Patch::initExchParticles(SmileiMPI* smpi, int ispec, PicParams& params, int tnum)
//{
//    Particles &cuParticles = (*vecSpecies[ispec]->particles);
//    
//    std::vector< std::vector<int> >* indexes_of_particles_to_exchange_per_thd = &vecSpecies[ispec]->indexes_of_particles_to_exchange_per_thd;
//    std::vector<int>*                indexes_of_particles_to_exchange         = &vecSpecies[ispec]->indexes_of_particles_to_exchange;
//    
//    /********************************************************************************/
//    // Build lists of indexes of particle to exchange per neighbor
//    // Computed from indexes_of_particles_to_exchange computed during particles' BC
//    /********************************************************************************/
//    (*indexes_of_particles_to_exchange).clear();
//        
//    int tmp = 0;
//    for (int tid=0 ; tid < (int)indexes_of_particles_to_exchange_per_thd->size() ; tid++)
//	tmp += ((*indexes_of_particles_to_exchange_per_thd)[tid]).size();
//    (*indexes_of_particles_to_exchange).resize( tmp );
//        
//    int k=0;
//    for (int tid=0 ; tid < (int)indexes_of_particles_to_exchange_per_thd->size() ; tid++) {
//	for (int ipart = 0 ; ipart < (int) ((*indexes_of_particles_to_exchange_per_thd)[tid]).size() ; ipart++ ) {
//	    (*indexes_of_particles_to_exchange)[k] =  (*indexes_of_particles_to_exchange_per_thd)[tid][ipart] ;
//	    k++;
//	}
//	((*indexes_of_particles_to_exchange_per_thd))[tid].clear();
//    }
//    sort( (*indexes_of_particles_to_exchange).begin(), (*indexes_of_particles_to_exchange).end() );
//        
//    int n_part_send = (*indexes_of_particles_to_exchange).size();
//    int n_part_recv;
//        
//    int iPart;
//    int n_particles;
//
//    // Define where particles are going 
//    for (int i=0 ; i<n_part_send ; i++) {
//	iPart = (*indexes_of_particles_to_exchange)[i];
//
//	if ( cuParticles.position(0,iPart) < min_local[0]) { 
//	    if ( cuParticles.position(1,iPart) < min_local[1]) {
//		vecSpecies[ispec]->specMPI.corner_buff_index_send[0][0].push_back( iPart );
//	    }
//	    else if ( cuParticles.position(1,iPart) >= max_local[1]) {
//		vecSpecies[ispec]->specMPI.corner_buff_index_send[0][1].push_back( iPart );
//	    }
//	    else {
//		vecSpecies[ispec]->specMPI.patch_buff_index_send[0][0].push_back( iPart );
//	    }
//	}
//	else if ( cuParticles.position(0,iPart) >= max_local[0]) { 
//	    if ( cuParticles.position(1,iPart) < min_local[1]) {
//		vecSpecies[ispec]->specMPI.corner_buff_index_send[1][0].push_back( iPart );
//	    }
//	    else if ( cuParticles.position(1,iPart) >= max_local[1]) {
//		vecSpecies[ispec]->specMPI.corner_buff_index_send[1][1].push_back( iPart );
//	    }
//	    else {
//		vecSpecies[ispec]->specMPI.patch_buff_index_send[0][1].push_back( iPart );
//	    }
//	}
//	else {
//	    if ( cuParticles.position(1,iPart) < min_local[1]) {
//		vecSpecies[ispec]->specMPI.patch_buff_index_send[1][0].push_back( iPart );
//	    }
//	    else if ( cuParticles.position(1,iPart) >= max_local[1]) {
//		vecSpecies[ispec]->specMPI.patch_buff_index_send[1][1].push_back( iPart );
//	    }
//	    else {
//		//If partciles is in but here, to be suppressed (= supp BC)
//	    }
//	}
//    }
//        
//    /********************************************************************************/
//    // Exchange number of particles to exchange to establish or not a communication
//    /********************************************************************************/
//    for (int iDim=0 ; iDim<2 ; iDim++) {
//	for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
//	    if (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {
//		//n_part_send = (vecSpecies[ispec]->specMPI.patch_buff_index_send[iDim][iNeighbor]).size();
//		vecSpecies[ispec]->specMPI.patch_buff_index_send_sz[iDim][iNeighbor] = (vecSpecies[ispec]->specMPI.patch_buff_index_send[iDim][iNeighbor]).size();
//		//int tag = buildtag( 1, hindex, neighbor_[iDim][iNeighbor]);
//		int tag = buildtag( hindex, iDim+1, iNeighbor+3 );
//		//MPI_Isend( &(vecSpecies[ispec]->specMPI.patch_buff_index_send_sz[iDim][iNeighbor]), 1, MPI_INT, 0, tag, MPI_COMM_SELF, &(vecSpecies[ispec]->specMPI.patch_srequest[iDim][iNeighbor]) );
//		MPI_Isend( &(vecSpecies[ispec]->specMPI.patch_buff_index_send_sz[iDim][iNeighbor]), 1, MPI_INT, MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &(vecSpecies[ispec]->specMPI.patch_srequest[iDim][iNeighbor]) );
//		//cout << hindex << " will sent " << n_part_send << " to " << neighbor_[iDim][iNeighbor] << " with tag " << tag << endl;
//	    } // END of Send
//	    else
//		n_part_send = 0;
//	    if (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
//		vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[iDim][(iNeighbor+1)%2] = 0;
//		//int tag = buildtag( 1, neighbor_[iDim][(iNeighbor+1)%2], hindex);
//		int tag = buildtag( neighbor_[iDim][(iNeighbor+1)%2], iDim+1, iNeighbor+3 );
//		//MPI_Irecv( &(vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[iDim][(iNeighbor+1)%2]), 1, MPI_INT, 0, tag, MPI_COMM_SELF, &(vecSpecies[ispec]->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]) );
//		MPI_Irecv( &(vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[iDim][(iNeighbor+1)%2]), 1, MPI_INT, MPI_neighbor_[iDim][(iNeighbor+1)%2], tag, MPI_COMM_WORLD, &(vecSpecies[ispec]->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]) );
//		//cout << hindex << " will recv from " << neighbor_[iDim][(iNeighbor+1)%2] << " with tag " << tag << endl;
//	    }
//	    else 
//		vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[iDim][(iNeighbor+1)%2] = 0;
//	}
//    }
//
//    for (int iDim=0 ; iDim<2 ; iDim++) {
//	for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
//	    if (corner_neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {
//		//n_part_send = (vecSpecies[ispec]->specMPI.corner_buff_index_send[iDim][iNeighbor]).size();
//		vecSpecies[ispec]->specMPI.corner_buff_index_send_sz[iDim][iNeighbor] = (vecSpecies[ispec]->specMPI.corner_buff_index_send[iDim][iNeighbor]).size();
//		//int tag = buildtag( 1, hindex, corner_neighbor_[iDim][iNeighbor]);
//		int tag = buildtag( 1, hindex, iDim+5, iNeighbor+7 );
//		//MPI_Isend( &(vecSpecies[ispec]->specMPI.corner_buff_index_send_sz[iDim][iNeighbor]), 1, MPI_INT, 0, tag, MPI_COMM_SELF, &(vecSpecies[ispec]->specMPI.corner_srequest[iDim][iNeighbor]) );
//		MPI_Isend( &(vecSpecies[ispec]->specMPI.corner_buff_index_send_sz[iDim][iNeighbor]), 1, MPI_INT, MPI_corner_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &(vecSpecies[ispec]->specMPI.corner_srequest[iDim][iNeighbor]) );
//		//cout << hindex << " will sent " << n_part_send << " to " << corner_neighbor_[iDim][iNeighbor] << " with tag " << tag << endl;
//	    } // END of Send
//	    else
//		n_part_send = 0;
//	    if (corner_neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
//		vecSpecies[ispec]->specMPI.corner_buff_index_recv_sz[iDim][(iNeighbor+1)%2] = 0;
//		//int tag = buildtag( 1, corner_neighbor_[iDim][(iNeighbor+1)%2], hindex);
//		int tag = buildtag( 1, corner_neighbor_[iDim][(iNeighbor+1)%2], (iDim+1)%2+5, iNeighbor+7 );
//		//MPI_Irecv( &(vecSpecies[ispec]->specMPI.corner_buff_index_recv_sz[iDim][(iNeighbor+1)%2]), 1, MPI_INT, 0, tag, MPI_COMM_SELF, &(vecSpecies[ispec]->specMPI.corner_rrequest[iDim][(iNeighbor+1)%2]) );
//		MPI_Irecv( &(vecSpecies[ispec]->specMPI.corner_buff_index_recv_sz[iDim][(iNeighbor+1)%2]), 1, MPI_INT, MPI_corner_neighbor_[iDim][(iNeighbor+1)%2], tag, MPI_COMM_WORLD, &(vecSpecies[ispec]->specMPI.corner_rrequest[iDim][(iNeighbor+1)%2]) );
//		//cout << hindex << " will recv from " << corner_neighbor_[iDim][(iNeighbor+1)%2] << " with tag " << tag << endl;
//	    }
//	    else 
//		vecSpecies[ispec]->specMPI.corner_buff_index_recv_sz[iDim][(iNeighbor+1)%2] = 0;
//	}
//    }
//
//} // initExchParticles


void Patch::initCommParticles(SmileiMPI* smpi, int ispec, PicParams& params, int tnum)
{
    Particles &cuParticles = (*vecSpecies[ispec]->particles);


    for (int iDim=0 ; iDim<2 ; iDim++) {
	for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
	    vecSpecies[ispec]->specMPI.patchVectorRecv[iDim][iNeighbor].initialize(0,2);
	    vecSpecies[ispec]->specMPI.patchVectorSend[iDim][iNeighbor].initialize(0,2);
	    vecSpecies[ispec]->specMPI.cornerVectorRecv[iDim][iNeighbor].initialize(0,2);
	    vecSpecies[ispec]->specMPI.cornerVectorSend[iDim][iNeighbor].initialize(0,2);
	}
    }


    /********************************************************************************/
    // Wait for end of communications over number of particles
    /********************************************************************************/
    for (int iDim=0 ; iDim<2 ; iDim++) {
	for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
	    MPI_Status sstat    [2];
	    MPI_Status rstat    [2];
	    if (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {
		MPI_Wait( &(vecSpecies[ispec]->specMPI.patch_srequest[iDim][iNeighbor]), &(sstat[iNeighbor]) );
		//cout << hindex << " has sent to " << neighbor_[iDim][iNeighbor]<< " using patch " << sstat[iNeighbor].MPI_TAG << endl;
	    }
	    if (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
		MPI_Wait( &(vecSpecies[ispec]->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]), &(rstat[(iNeighbor+1)%2]) );
		//cout << hindex << " will recv " << vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[iDim][(iNeighbor+1)%2] << " particles from " << neighbor_[iDim][(iNeighbor+1)%2]  << " with tag " << rstat[(iNeighbor+1)%2].MPI_TAG << endl;
		if (vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[iDim][(iNeighbor+1)%2]!=0) {
		    vecSpecies[ispec]->specMPI.patchVectorRecv[iDim][(iNeighbor+1)%2].initialize( vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[iDim][(iNeighbor+1)%2], cuParticles.dimension());
		}
	    }
	}
    }
    for (int iDim=0 ; iDim<2 ; iDim++) {
	for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
	    MPI_Status sstat    [2];
	    MPI_Status rstat    [2];
	    if (corner_neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {
		MPI_Wait( &(vecSpecies[ispec]->specMPI.corner_srequest[iDim][iNeighbor]), &(sstat[iNeighbor]) );
	    }
	    if (corner_neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
		MPI_Wait( &(vecSpecies[ispec]->specMPI.corner_rrequest[iDim][(iNeighbor+1)%2]), &(rstat[(iNeighbor+1)%2]) );
		//cout << hindex << " will recv " << vecSpecies[ispec]->specMPI.corner_buff_index_recv_sz[iDim][(iNeighbor+1)%2] << " particles from " << corner_neighbor_[iDim][(iNeighbor+1)%2] << " with tag " << rstat[(iNeighbor+1)%2].MPI_TAG << endl;
		if (vecSpecies[ispec]->specMPI.corner_buff_index_recv_sz[iDim][(iNeighbor+1)%2]!=0) {
		    vecSpecies[ispec]->specMPI.cornerVectorRecv[iDim][(iNeighbor+1)%2].initialize( vecSpecies[ispec]->specMPI.corner_buff_index_recv_sz[iDim][(iNeighbor+1)%2], cuParticles.dimension());
		}
	    }
	}
    }
    
    /********************************************************************************/
    // Proceed to effective Particles' communications
    /********************************************************************************/

    // Number of properties per particles = nDim_Particles + 3 + 1 + 1
    int nbrOfProp( 7 );
    MPI_Datatype typePartSend, typePartRecv;

    int n_part_send, n_part_recv;

    for (int iDim=0 ; iDim<2 ; iDim++) {
	for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
                
	    // n_part_send : number of particles to send to current neighbor
	    n_part_send = (vecSpecies[ispec]->specMPI.patch_buff_index_send[iDim][iNeighbor]).size();
	    if ( (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) && (n_part_send!=0) ) {
		double x_max = params.cell_length[iDim]*( params.n_space_global[iDim] );
		for (int iPart=0 ; iPart<n_part_send ; iPart++) {
		    if (smpi->periods_[iDim]==1) {
			// Enabled periodicity
			if ( ( iNeighbor==0 ) &&  (Pcoordinates[iDim] == 0 ) &&( cuParticles.position(iDim,vecSpecies[ispec]->specMPI.patch_buff_index_send[iDim][iNeighbor][iPart]) < 0. ) ) {
			    cuParticles.position(iDim,vecSpecies[ispec]->specMPI.patch_buff_index_send[iDim][iNeighbor][iPart])     += x_max;
			}
			else if ( ( iNeighbor==1 ) &&  (Pcoordinates[iDim] == params.number_of_patches[iDim]-1 ) && ( cuParticles.position(iDim,vecSpecies[ispec]->specMPI.patch_buff_index_send[iDim][iNeighbor][iPart]) >= x_max ) ) {
			    cuParticles.position(iDim,vecSpecies[ispec]->specMPI.patch_buff_index_send[iDim][iNeighbor][iPart])     -= x_max;
			}
		    }
		    cuParticles.cp_particle(vecSpecies[ispec]->specMPI.patch_buff_index_send[iDim][iNeighbor][iPart], vecSpecies[ispec]->specMPI.patchVectorSend[iDim][iNeighbor]);
		}

		//int tag = buildtag( 2, hindex, neighbor_[iDim][iNeighbor]);
		int tag = buildtag( hindex, iDim+1, iNeighbor+3 );
		typePartSend = smpi->createMPIparticles( &(vecSpecies[ispec]->specMPI.patchVectorSend[iDim][iNeighbor]), nbrOfProp );
		//MPI_Isend( &((vecSpecies[ispec]->specMPI.patchVectorSend[iDim][iNeighbor]).position(0,0)), 1, typePartSend, 0, tag, MPI_COMM_SELF, &(vecSpecies[ispec]->specMPI.patch_srequest[iDim][iNeighbor]) );
		MPI_Isend( &((vecSpecies[ispec]->specMPI.patchVectorSend[iDim][iNeighbor]).position(0,0)), 1, typePartSend, MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &(vecSpecies[ispec]->specMPI.patch_srequest[iDim][iNeighbor]) );
		//cout << hindex << " really send " << n_part_send << " to " << neighbor_[iDim][iNeighbor] << " with tag " << tag << endl;
		MPI_Type_free( &typePartSend );

	    } // END of Send
                
	    n_part_recv = vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[iDim][(iNeighbor+1)%2];
	    if ( (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
		typePartRecv = smpi->createMPIparticles( &(vecSpecies[ispec]->specMPI.patchVectorRecv[iDim][(iNeighbor+1)%2]), nbrOfProp );
		int tag = buildtag( neighbor_[iDim][(iNeighbor+1)%2], iDim+1, iNeighbor+3 );
		//int tag = buildtag( 2, neighbor_[iDim][(iNeighbor+1)%2], hindex);
		//MPI_Irecv( &((vecSpecies[ispec]->specMPI.patchVectorRecv[iDim][(iNeighbor+1)%2]).position(0,0)), 1, typePartRecv,  0, tag, MPI_COMM_SELF, &(vecSpecies[ispec]->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]) );
		MPI_Irecv( &((vecSpecies[ispec]->specMPI.patchVectorRecv[iDim][(iNeighbor+1)%2]).position(0,0)), 1, typePartRecv, MPI_neighbor_[iDim][(iNeighbor+1)%2], tag, MPI_COMM_WORLD, &(vecSpecies[ispec]->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]) );
		//cout << hindex << " will really recv " << n_part_recv << " from " << neighbor_[iDim][(iNeighbor+1)%2] << " with tag " << tag << endl;
		MPI_Type_free( &typePartRecv );

	    } // END of Recv
                
	} // END for iNeighbor
    } // END for iDim


    for (int iDim=0 ; iDim<2 ; iDim++) {
	for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
                
	    // n_part_send : number of particles to send to current neighbor
	    n_part_send = (vecSpecies[ispec]->specMPI.corner_buff_index_send[iDim][iNeighbor]).size();
	    if ( (corner_neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) && (n_part_send!=0) ) {
		double x_max = params.cell_length[(iDim+1)%2]*( params.n_space_global[(iDim+1)%2] );
		double y_max = params.cell_length[iDim]*( params.n_space_global[iDim] );
		for (int iPart=0 ; iPart<n_part_send ; iPart++) {
		    if (smpi->periods_[iDim]==1) {
			// Enabled periodicity
			if ( (Pcoordinates[(iDim+1)%2] == 0 ) &&( cuParticles.position((iDim+1)%2,vecSpecies[ispec]->specMPI.corner_buff_index_send[iDim][iNeighbor][iPart]) < 0. ) ) {
			    cuParticles.position((iDim+1)%2,vecSpecies[ispec]->specMPI.corner_buff_index_send[iDim][iNeighbor][iPart])     += x_max;
			}
			if ( (Pcoordinates[iDim] == 0 ) &&( cuParticles.position(iDim,vecSpecies[ispec]->specMPI.corner_buff_index_send[iDim][iNeighbor][iPart]) < 0. ) ) {
			    cuParticles.position(iDim,vecSpecies[ispec]->specMPI.corner_buff_index_send[iDim][iNeighbor][iPart])     += y_max;
			}
			if ( (Pcoordinates[(iDim+1)%2] == params.number_of_patches[(iDim+1)%2]-1 ) && ( cuParticles.position((iDim+1)%2,vecSpecies[ispec]->specMPI.corner_buff_index_send[iDim][iNeighbor][iPart]) >= x_max ) ) {
			    cuParticles.position((iDim+1)%2,vecSpecies[ispec]->specMPI.corner_buff_index_send[iDim][iNeighbor][iPart])     -= x_max;
			}
			if ( (Pcoordinates[iDim] == params.number_of_patches[iDim]-1 ) && ( cuParticles.position(iDim,vecSpecies[ispec]->specMPI.corner_buff_index_send[iDim][iNeighbor][iPart]) >= y_max ) ) {
			    cuParticles.position(iDim,vecSpecies[ispec]->specMPI.corner_buff_index_send[iDim][iNeighbor][iPart])     -= y_max;
			}

		    }
		    cuParticles.cp_particle(vecSpecies[ispec]->specMPI.corner_buff_index_send[iDim][iNeighbor][iPart], vecSpecies[ispec]->specMPI.cornerVectorSend[iDim][iNeighbor]);
		}

		typePartSend = smpi->createMPIparticles( &(vecSpecies[ispec]->specMPI.cornerVectorSend[iDim][iNeighbor]), nbrOfProp );
		//int tag = buildtag( 2, hindex, corner_neighbor_[iDim][iNeighbor]);
		int tag = buildtag( 1, hindex, iDim, iNeighbor );
		//MPI_Isend( &((vecSpecies[ispec]->specMPI.cornerVectorSend[iDim][iNeighbor]).position(0,0)), 1, typePartSend, 0, tag, MPI_COMM_SELF, &(vecSpecies[ispec]->specMPI.corner_srequest[iDim][iNeighbor]) );
		MPI_Isend( &((vecSpecies[ispec]->specMPI.cornerVectorSend[iDim][iNeighbor]).position(0,0)), 1, typePartSend, MPI_corner_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &(vecSpecies[ispec]->specMPI.corner_srequest[iDim][iNeighbor]) );
		//cout << hindex << " really send " << n_part_send << " to " << corner_neighbor_[iDim][iNeighbor] << " with tag " << tag << endl;
		MPI_Type_free( &typePartSend );

	    } // END of Send
                
	    n_part_recv = vecSpecies[ispec]->specMPI.corner_buff_index_recv_sz[iDim][(iNeighbor+1)%2];
	    if ( (corner_neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
		typePartRecv = smpi->createMPIparticles( &(vecSpecies[ispec]->specMPI.cornerVectorRecv[iDim][(iNeighbor+1)%2]), nbrOfProp );
		//int tag = buildtag( 2, corner_neighbor_[iDim][(iNeighbor+1)%2], hindex);
		int tag = buildtag( 1, corner_neighbor_[iDim][(iNeighbor+1)%2], (iDim+1)%2, iNeighbor );
		//MPI_Irecv( &((vecSpecies[ispec]->specMPI.cornerVectorRecv[iDim][(iNeighbor+1)%2]).position(0,0)), 1, typePartRecv, 0, tag, MPI_COMM_SELF, &(vecSpecies[ispec]->specMPI.corner_rrequest[iDim][(iNeighbor+1)%2]) );
		MPI_Irecv( &((vecSpecies[ispec]->specMPI.cornerVectorRecv[iDim][(iNeighbor+1)%2]).position(0,0)), 1, typePartRecv, MPI_corner_neighbor_[iDim][(iNeighbor+1)%2], tag, MPI_COMM_WORLD, &(vecSpecies[ispec]->specMPI.corner_rrequest[iDim][(iNeighbor+1)%2]) );
		//cout << hindex << " will really recv from " << corner_neighbor_[iDim][(iNeighbor+1)%2] << " with tag " << tag << endl;
		MPI_Type_free( &typePartRecv );

	    } // END of Recv
                
	} // END for iNeighbor
    } // END for iDim


} // initCommParticles

void Patch::finalizeCommParticles(SmileiMPI* smpi, int ispec, PicParams& params, int tnum)
{
    Particles &cuParticles = (*vecSpecies[ispec]->particles);


    //std::vector< std::vector<int> >* indexes_of_particles_to_exchange_per_thd = &vecSpecies[ispec]->indexes_of_particles_to_exchange_per_thd;
    std::vector<int>*                indexes_of_particles_to_exchange         = &vecSpecies[ispec]->indexes_of_particles_to_exchange;

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
    for (int iDim=0 ; iDim<2 ; iDim++) {
	for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
	    MPI_Status sstat    [2];
	    MPI_Status rstat    [2];
                
	    n_part_send = vecSpecies[ispec]->specMPI.patch_buff_index_send[iDim][iNeighbor].size();
	    n_part_recv = vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[iDim][(iNeighbor+1)%2];
                
	    if ( (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) && (n_part_send!=0) ) {
		MPI_Wait( &(vecSpecies[ispec]->specMPI.patch_srequest[iDim][iNeighbor]), &(sstat[iNeighbor]) );
		// clean vecSpecies[ispec]->specMPI.patchVectorSend
		//vecSpecies[ispec]->specMPI.patchVectorSend[iDim][iNeighbor].erase_particle_trail(0);
	    }
	    if ( (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
		MPI_Wait( &(vecSpecies[ispec]->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]), &(rstat[(iNeighbor+1)%2]) );     
		//cout << hindex << " recv from " << neighbor_[iDim][(iNeighbor+1)%2] << endl;
	    }
	}
    }

    for (int iDim=0 ; iDim<2 ; iDim++) {
	for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
                
	    MPI_Status sstat    [2];
	    MPI_Status rstat    [2];

	    n_part_send = vecSpecies[ispec]->specMPI.corner_buff_index_send[iDim][iNeighbor].size();
	    n_part_recv = vecSpecies[ispec]->specMPI.corner_buff_index_recv_sz[iDim][(iNeighbor+1)%2];
                
	    if ( (corner_neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) && (n_part_send!=0) ) {
		MPI_Wait( &(vecSpecies[ispec]->specMPI.corner_srequest[iDim][iNeighbor]), &(sstat[iNeighbor]) );
		// clean cornerVectorSend
		//vecSpecies[ispec]->specMPI.cornerVectorSend[iDim][iNeighbor].erase_particle_trail(0);
	    }
	    if ( (corner_neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
		MPI_Wait( &(vecSpecies[ispec]->specMPI.corner_rrequest[iDim][(iNeighbor+1)%2]), &(rstat[(iNeighbor+1)%2]) );
		//cout << hindex << " recv from " << corner_neighbor_[iDim][(iNeighbor+1)%2] << endl;
	    }
	}
    }

    /********************************************************************************/
    // Clean lists of indexes of particle to exchange per neighbor
    /********************************************************************************/
    for (int iDim=0 ; iDim<2 ; iDim++) {
	for (int i=0 ; i<nbNeighbors_ ; i++) {
	    vecSpecies[ispec]->specMPI.patch_buff_index_send[iDim][i].clear();
	    vecSpecies[ispec]->specMPI.corner_buff_index_send[iDim][i].clear();
	}
    }
            
        
    /********************************************************************************/
    // Delete Particles included in buff_send/buff_recv
    /********************************************************************************/
    cleanup_sent_particles(ispec, indexes_of_particles_to_exchange);


    // Delete useless Particles
    //Theoretically, not even necessary to do anything as long you use bmax as the end of your iterator on particles.
    //Nevertheless, you might want to free memory and have the actual number of particles
    //really equal to the size of the vector. So we do:
    cuParticles.erase_particle_trail((*cubmax).back());
    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



    for (int iDim=0 ; iDim<2 ; iDim++) {
	for (unsigned int j=0; j<(*cubmax).size()+1 ;j++){
	    shift[j]=0;
	}

	//********************************************************************************/
	// Copy newly arrived particles back to the vector
	// WARNING: very different behaviour depending on which dimension particles are coming from.
	/********************************************************************************/
	//We first evaluate how many particles arrive in each bin.
	if (iDim==1) {
	    //1) Count particles coming from south and north
	    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
		n_part_recv = vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[iDim][iNeighbor];
		for (unsigned int j=0; j<n_part_recv ;j++){
		    ii = int((vecSpecies[ispec]->specMPI.patchVectorRecv[iDim][iNeighbor].position(0,j)-min_local[0])/dbin);//bin in which the particle goes.
		    shift[ii+1]++; // It makes the next bins shift.
		}
	    }
	}
	if (iDim==0) {
	    //2) Add particles coming from west and east
	    shift[1] += vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[iDim][0];//Particles coming from south all go to bin 0 and shift all the other bins.
	    shift[(*cubmax).size()] += vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[iDim][1];//Used only to count the total number of particles arrived.
	}


	//Evaluation of the necessary shift of all bins.
	//Must be done sequentially
	for (unsigned int j=1; j<(*cubmax).size()+1;j++){ //bin 0 is not shifted.Last element of shift stores total number of arriving particles.
	    shift[j]+=shift[j-1];
	}
	//Make room for new particles
	cuParticles.initialize( cuParticles.size()+shift[(*cubmax).size()], cuParticles.dimension() );
        
	//Shift bins, must be done sequentially
	for (unsigned int j=(*cubmax).size()-1; j>=1; j--){
	    n_particles = (*cubmax)[j]-(*cubmin)[j]; //Nbr of particle in this bin
	    nmove = min(n_particles,shift[j]); //Nbr of particles to move
	    lmove = max(n_particles,shift[j]); //How far particles must be shifted
	    if (nmove>0) cuParticles.overwrite_part2D((*cubmin)[j], (*cubmin)[j]+lmove, nmove);
	    (*cubmin)[j] += shift[j];
	    (*cubmax)[j] += shift[j];
	}
        
	//Space has been made now to write the arriving particles into the correct bins
	//iDim == 0  is the easy case, when particles arrive either in first or last bin.
	if (iDim==0) {
	    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
		n_part_recv = vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[iDim][iNeighbor];
		if ( (neighbor_[0][iNeighbor]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
		    ii = iNeighbor*((*cubmax).size()-1);//0 if iNeighbor=0(particles coming from West) and (*cubmax).size()-1 otherwise.
		    vecSpecies[ispec]->specMPI.patchVectorRecv[iDim][iNeighbor].overwrite_part2D(0, cuParticles,(*cubmax)[ii],n_part_recv);
		    (*cubmax)[ii] += n_part_recv ;
		}
	    }
	}
	//iDim == 1) is the difficult case, when particles can arrive in any bin.
	if (iDim==1) {
	    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
		n_part_recv = vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[iDim][iNeighbor];
		if ( (neighbor_[1][iNeighbor]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
		    for(unsigned int j=0; j<n_part_recv; j++){
			ii = int((vecSpecies[ispec]->specMPI.patchVectorRecv[iDim][iNeighbor].position(0,j)-min_local[0])/dbin);//bin in which the particle goes.
			vecSpecies[ispec]->specMPI.patchVectorRecv[iDim][iNeighbor].overwrite_part2D(j, cuParticles,(*cubmax)[ii]);
			(*cubmax)[ii] ++ ;
		    }
		}
	    }
	}
    } // End for iDim

    // ------------------ CORNERS ------------------
    for (int iDim=0 ; iDim<2 ; iDim++) {
	for (unsigned int j=0; j<(*cubmax).size()+1 ;j++){
	    shift[j]=0;
	}

	//********************************************************************************/
	// Copy newly arrived particles back to the vector
	// WARNING: very different behaviour depending on which dimension particles are coming from.
	/********************************************************************************/
	//We first evaluate how many particles arrive in each bin.
	//2) Add particles coming from west and east
	if (iDim==0) {
	    shift[1] += vecSpecies[ispec]->specMPI.corner_buff_index_recv_sz[iDim][0];//Particles coming from south all go to bin 0 and shift all the other bins.
	    shift[1] += vecSpecies[ispec]->specMPI.corner_buff_index_recv_sz[iDim][1];//Particles coming from south all go to bin 0 and shift all the other bins.
	}
	if (iDim==1) {
	    shift[(*cubmax).size()] += vecSpecies[ispec]->specMPI.corner_buff_index_recv_sz[iDim][0];//Used only to count the total number of particles arrived.
	    shift[(*cubmax).size()] += vecSpecies[ispec]->specMPI.corner_buff_index_recv_sz[iDim][1];//Used only to count the total number of particles arrived.
	}

	//Evaluation of the necessary shift of all bins.
	//Must be done sequentially
	for (unsigned int j=1; j<(*cubmax).size()+1;j++){ //bin 0 is not shifted.Last element of shift stores total number of arriving particles.
	    shift[j]+=shift[j-1];
	}
	//Make room for new particles
	cuParticles.initialize( cuParticles.size()+shift[(*cubmax).size()], cuParticles.dimension() );
        
	//Shift bins, must be done sequentially
	for (unsigned int j=(*cubmax).size()-1; j>=1; j--){
	    n_particles = (*cubmax)[j]-(*cubmin)[j]; //Nbr of particle in this bin
	    nmove = min(n_particles,shift[j]); //Nbr of particles to move
	    lmove = max(n_particles,shift[j]); //How far particles must be shifted
	    if (nmove>0) cuParticles.overwrite_part2D((*cubmin)[j], (*cubmin)[j]+lmove, nmove);
	    (*cubmin)[j] += shift[j];
	    (*cubmax)[j] += shift[j];
	}
        
	//Space has been made now to write the arriving particles into the correct bins
	//Corner particles arrive either in first or last bin.
	for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
	    n_part_recv = vecSpecies[ispec]->specMPI.corner_buff_index_recv_sz[iDim][iNeighbor];
	    if ( (corner_neighbor_[0][iNeighbor]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
		ii = iDim*((*cubmax).size()-1);//0 if iDim=0(particles coming from West) and (*cubmax).size()-1 otherwise.
		vecSpecies[ispec]->specMPI.cornerVectorRecv[iDim][iNeighbor].overwrite_part2D(0, cuParticles,(*cubmax)[ii],n_part_recv);
		(*cubmax)[ii] += n_part_recv ;
	    }
	}
    } // End for iDim
} // finalizeCommParticles


void Patch::initExchParticles(SmileiMPI* smpi, int ispec, PicParams& params, int tnum, VectorPatch * vecPatch)
{
    Particles &cuParticles = (*vecSpecies[ispec]->particles);
    int ndim = params.nDim_field;
    int idim,check;
    std::vector<int>*                indexes_of_particles_to_exchange         = &vecSpecies[ispec]->indexes_of_particles_to_exchange;
    double xmax[3]; 
    
    for (int iDim=0 ; iDim < ndim ; iDim++){
        xmax[iDim] = params.cell_length[iDim]*( params.n_space_global[iDim] );
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
            vecSpecies[ispec]->specMPI.patchVectorRecv[iDim][iNeighbor].initialize(0,2);
            vecSpecies[ispec]->specMPI.patchVectorSend[iDim][iNeighbor].initialize(0,2);
	    vecSpecies[ispec]->specMPI.patch_buff_index_send[idim][iNeighbor].clear();
        }
    }

    

     sort( (*indexes_of_particles_to_exchange).begin(), (*indexes_of_particles_to_exchange).end() );
 
 
    int n_part_send = (*indexes_of_particles_to_exchange).size();
    //int n_part_recv;
        
    int iPart;
    int n_particles;

    // Define where particles are going 
    //diagonalParticles.initialize(0,cuParticles.dimension());
    //Put particles in the send buffer it belongs to. Priority to lower dimensions.
    for (int i=0 ; i<n_part_send ; i++) {
	iPart = (*indexes_of_particles_to_exchange)[i];
        check = 0;
        idim = 0;

        //Put indexes of particles in the first direction they will be exchanged and correct their position according to periodicity for the first exchange only.
        while (check == 0 && idim<ndim){
	    if ( cuParticles.position(idim,iPart) < min_local[idim]) { 
	        vecSpecies[ispec]->specMPI.patch_buff_index_send[idim][0].push_back( iPart );
		if (smpi->periods_[idim]==1 && Pcoordinates[idim] == 0) {
	            cuParticles.position(idim,iPart)     += xmax[idim];
                }
                check = 1;
                cout << "sening part along dim -" << idim << endl;
	    }
	    else if ( cuParticles.position(idim,iPart) >= max_local[idim]) { 
	        vecSpecies[ispec]->specMPI.patch_buff_index_send[idim][1].push_back( iPart );
		if (smpi->periods_[idim]==1 && Pcoordinates[idim] == params.number_of_patches[idim]-1) {
	            cuParticles.position(idim,iPart)     -= xmax[idim];
                }
                check = 1;
                cout << "sening part along dim +" << idim << endl;
	    }
            idim++;
        }
    }
    //cout << hindex << " has " << n_part_send << " to send" << endl;
     

} // initExchParticles(... iDim)

void Patch::initCommParticles(SmileiMPI* smpi, int ispec, PicParams& params, int tnum, int iDim, VectorPatch * vecPatch)
{

    Particles &cuParticles = (*vecSpecies[ispec]->particles);
    int nbrOfProp( 7 );
    MPI_Datatype typePartSend, typePartRecv;

    int n_part_send, n_part_recv;

    int h0 = (*vecPatch)(0)->hindex;
    /********************************************************************************/
    // Exchange number of particles to exchange to establish or not a communication
    /********************************************************************************/
    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
	if (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {
	    vecSpecies[ispec]->specMPI.patch_buff_index_send_sz[iDim][iNeighbor] = (vecSpecies[ispec]->specMPI.patch_buff_index_send[iDim][iNeighbor]).size();
	    if (is_a_MPI_neighbor(iDim, iNeighbor)) {
		//cout << " send through MPI" << endl;
		//int tag = buildtag( 1, hindex, neighbor_[iDim][iNeighbor]);
		int tag = buildtag( hindex, iDim+1, iNeighbor+3 );
                if (vecSpecies[ispec]->specMPI.patch_buff_index_send_sz[iDim][iNeighbor]!=0){
                    cout << smpi->smilei_rk << "send " << vecSpecies[ispec]->specMPI.patch_buff_index_send_sz[iDim][iNeighbor] << " to " << MPI_neighbor_[iDim][iNeighbor] << "with tag " << tag << endl;
                }
		MPI_Isend( &(vecSpecies[ispec]->specMPI.patch_buff_index_send_sz[iDim][iNeighbor]), 1, MPI_INT, MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &(vecSpecies[ispec]->specMPI.patch_srequest[iDim][iNeighbor]) );}
	    else {
		//cout << hindex << " send " << neighbor_[iDim][iNeighbor]- h0 << " with " << vecSpecies[ispec]->specMPI.patch_buff_index_send_sz[iDim][iNeighbor] << " particles" << endl;
		(*vecPatch)( neighbor_[iDim][iNeighbor]- h0 )->vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[iDim][(iNeighbor+1)%2] = vecSpecies[ispec]->specMPI.patch_buff_index_send_sz[iDim][iNeighbor];
		//cout << neighbor_[iDim][iNeighbor]- h0 << " initialize with " << vecSpecies[ispec]->specMPI.patch_buff_index_send_sz[iDim][iNeighbor] << " particles" << endl;
		//(*vecPatch)( neighbor_[iDim][iNeighbor]- h0 )->vecSpecies[ispec]->specMPI.patchVectorRecv[iDim][(iNeighbor+1)%2].initialize( vecSpecies[ispec]->specMPI.patch_buff_index_send_sz[iDim][iNeighbor], cuParticles.dimension());

	    }
	} // END of Send

	if (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
	    if (is_a_MPI_neighbor(iDim, (iNeighbor+1)%2)) {
		//cout << " recv through MPI" << endl;
		//Is this init to 0 really necessary ??
		vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[iDim][(iNeighbor+1)%2] = 0;
		//int tag = buildtag( 1, neighbor_[iDim][(iNeighbor+1)%2], hindex);
		int tag = buildtag( neighbor_[iDim][(iNeighbor+1)%2], iDim+1, iNeighbor+3 );
                cout << smpi->smilei_rk << "receives from " << MPI_neighbor_[iDim][(iNeighbor+1)%2] << " with tag " << tag << endl;
		MPI_Irecv( &(vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[iDim][(iNeighbor+1)%2]), 1, MPI_INT, MPI_neighbor_[iDim][(iNeighbor+1)%2], tag, MPI_COMM_WORLD, &(vecSpecies[ispec]->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]) );
	    }
	}
	else
	    vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[iDim][(iNeighbor+1)%2] = 0;

    }

} // initCommParticles(... iDim)

void Patch::finalizeCommParticles(SmileiMPI* smpi, int ispec, PicParams& params, int tnum, int iDim, VectorPatch * vecPatch)
{
    int ndim = params.nDim_field;
    int idim, check;
    int nbrOfProp( 7 );
    MPI_Datatype typePartSend, typePartRecv;

    Particles &cuParticles = (*vecSpecies[ispec]->particles);
    double xmax[3]; 
    for (int idim=0 ; idim < ndim ; idim++)
        xmax[iDim] = params.cell_length[iDim]*( params.n_space_global[iDim] );


    std::vector<int>*                indexes_of_particles_to_exchange         = &vecSpecies[ispec]->indexes_of_particles_to_exchange;

    std::vector<int>* cubmin = &vecSpecies[ispec]->bmin;
    std::vector<int>* cubmax = &vecSpecies[ispec]->bmax;

    int nmove,lmove,ii; // local, OK
    int shift[(*cubmax).size()+1];//how much we need to shift each bin in order to leave room for the new particles
    double dbin;
        
    dbin = params.cell_length[0]*params.clrw; //width of a bin.
    int h0 = (*vecPatch)(0)->hindex;
    for (unsigned int j=0; j<(*cubmax).size()+1 ;j++){
      shift[j]=0;
    }

    int n_part_send, n_part_recv, n_particles;

    /********************************************************************************/
    // Wait for end of communications over number of particles
    /********************************************************************************/
    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
	MPI_Status sstat    [2];
	MPI_Status rstat    [2];
	if (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {
	    if (is_a_MPI_neighbor(iDim, iNeighbor))
		MPI_Wait( &(vecSpecies[ispec]->specMPI.patch_srequest[iDim][iNeighbor]), &(sstat[iNeighbor]) );
	}
	if (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
	    if (is_a_MPI_neighbor(iDim, (iNeighbor+1)%2))  {
		MPI_Wait( &(vecSpecies[ispec]->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]), &(rstat[(iNeighbor+1)%2]) );
		if (vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[iDim][(iNeighbor+1)%2]!=0) {
		    //cout << hindex << " initialize with " << vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[iDim][(iNeighbor+1)%2] << " particles" << endl;
		    vecSpecies[ispec]->specMPI.patchVectorRecv[iDim][(iNeighbor+1)%2].initialize( vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[iDim][(iNeighbor+1)%2], cuParticles.dimension());
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
	n_part_send = (vecSpecies[ispec]->specMPI.patch_buff_index_send[iDim][iNeighbor]).size();
	if ( (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) && (n_part_send!=0) ) {
            cout << "npart s is non zero !!! dim= " <<iDim << " nparts= "<<n_part_send << endl;
	    double x_max = params.cell_length[iDim]*( params.n_space_global[iDim] );
	    for (int iPart=0 ; iPart<n_part_send ; iPart++) {
		if (smpi->periods_[iDim]==1) {
		    // Enabled periodicity
		    if ( ( iNeighbor==0 ) &&  (Pcoordinates[iDim] == 0 ) &&( cuParticles.position(iDim,vecSpecies[ispec]->specMPI.patch_buff_index_send[iDim][iNeighbor][iPart]) < 0. ) ) {
			cuParticles.position(iDim,vecSpecies[ispec]->specMPI.patch_buff_index_send[iDim][iNeighbor][iPart])     += x_max;
		    }
		    else if ( ( iNeighbor==1 ) &&  (Pcoordinates[iDim] == params.number_of_patches[iDim]-1 ) && ( cuParticles.position(iDim,vecSpecies[ispec]->specMPI.patch_buff_index_send[iDim][iNeighbor][iPart]) >= x_max ) ) {
			cuParticles.position(iDim,vecSpecies[ispec]->specMPI.patch_buff_index_send[iDim][iNeighbor][iPart])     -= x_max;
		    }
		}
	    }
	    if (is_a_MPI_neighbor(iDim, iNeighbor)) {
	        for (int iPart=0 ; iPart<n_part_send ; iPart++) 
		    cuParticles.cp_particle(vecSpecies[ispec]->specMPI.patch_buff_index_send[iDim][iNeighbor][iPart], vecSpecies[ispec]->specMPI.patchVectorSend[iDim][iNeighbor]);
		//int tag = buildtag( 2, hindex, neighbor_[iDim][iNeighbor]);
		int tag = buildtag( hindex, iDim+1, iNeighbor+3 );
		typePartSend = smpi->createMPIparticles( &(vecSpecies[ispec]->specMPI.patchVectorSend[iDim][iNeighbor]), nbrOfProp );
		MPI_Isend( &((vecSpecies[ispec]->specMPI.patchVectorSend[iDim][iNeighbor]).position(0,0)), 1, typePartSend, MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &(vecSpecies[ispec]->specMPI.patch_srequest[iDim][iNeighbor]) );
		MPI_Type_free( &typePartSend );
	    }
	    else {
	        for (int iPart=0 ; iPart<n_part_send ; iPart++) 
		    cuParticles.cp_particle( vecSpecies[ispec]->specMPI.patch_buff_index_send[iDim][iNeighbor][iPart],((*vecPatch)( neighbor_[iDim][iNeighbor]- h0 )->vecSpecies[ispec]->specMPI.patchVectorRecv[iDim][(iNeighbor+1)%2]) );
		//(vecSpecies[ispec]->specMPI.patchVectorSend[iDim][iNeighbor]).cp_particles(0, n_part_send, ((*vecPatch)( neighbor_[iDim][iNeighbor]- h0 )->vecSpecies[ispec]->specMPI.patchVectorRecv[iDim][(iNeighbor+1)%2]), 0);

	    }

	} // END of Send
                
	n_part_recv = vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[iDim][(iNeighbor+1)%2];
	if ( (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
	    if (is_a_MPI_neighbor(iDim, (iNeighbor+1)%2)) {
		typePartRecv = smpi->createMPIparticles( &(vecSpecies[ispec]->specMPI.patchVectorRecv[iDim][(iNeighbor+1)%2]), nbrOfProp );
		//int tag = buildtag( 2, neighbor_[iDim][(iNeighbor+1)%2], hindex);
		int tag = buildtag( neighbor_[iDim][(iNeighbor+1)%2], iDim+1 ,iNeighbor+3 );
		MPI_Irecv( &((vecSpecies[ispec]->specMPI.patchVectorRecv[iDim][(iNeighbor+1)%2]).position(0,0)), 1, typePartRecv, MPI_neighbor_[iDim][(iNeighbor+1)%2], tag, MPI_COMM_WORLD, &(vecSpecies[ispec]->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]) );
		MPI_Type_free( &typePartRecv );
	    }

	} // END of Recv
                
    } // END for iNeighbor




    /********************************************************************************/
    // Wait for end of communications over Particles
    /********************************************************************************/
    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
	MPI_Status sstat    [2];
	MPI_Status rstat    [2];
                
	n_part_send = vecSpecies[ispec]->specMPI.patch_buff_index_send[iDim][iNeighbor].size();
	n_part_recv = vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[iDim][(iNeighbor+1)%2];
               
        if (iDim == ndim-1){
            //We have stored in indexes_of_particles_to_exchange the list of all particles that needs to be removed.
            cleanup_sent_particles(ispec, indexes_of_particles_to_exchange);
            (*indexes_of_particles_to_exchange).clear();
            cuParticles.erase_particle_trail((*cubmax).back());
        }

 
	if ( (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) && (n_part_send!=0) ) {
	    if (is_a_MPI_neighbor(iDim, iNeighbor))
		MPI_Wait( &(vecSpecies[ispec]->specMPI.patch_srequest[iDim][iNeighbor]), &(sstat[iNeighbor]) );
	}
	if ( (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
	    if (is_a_MPI_neighbor(iDim, (iNeighbor+1)%2))
		MPI_Wait( &(vecSpecies[ispec]->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]), &(rstat[(iNeighbor+1)%2]) );     

	    // Treat diagonalParticles
	    if (iDim < ndim-1){ // No need to treat diag particles at last dimension.
                //for (int iPart=0 ; iPart<n_part_recv; iPart++ ){
	        for (int iPart=n_part_recv-1 ; iPart>=0; iPart-- ) {
                    check = 0;
                    idim = iDim+1;//We check next dimension
                    while (check == 0 && idim<ndim){
                        //If particle not in the domain...
	                if ( (vecSpecies[ispec]->specMPI.patchVectorRecv[iDim][(iNeighbor+1)%2]).position(idim,iPart) < min_local[idim]) { 
                            //...Deal with periodicity...
	            	    if (smpi->periods_[idim]==1 && Pcoordinates[idim] == 0) {
	                        (vecSpecies[ispec]->specMPI.patchVectorRecv[iDim][(iNeighbor+1)%2]).position(idim,iPart)     += xmax[idim];
                            }
                            //... copy it at the back of the local particle vector ...
                            (vecSpecies[ispec]->specMPI.patchVectorRecv[iDim][(iNeighbor+1)%2]).cp_particle(iPart, cuParticles);
                            //...adjust bmax ...
                            (*cubmax)[(*cubmax).size()-1]++;
                            //... and add its index to the particles to be sent later...
	                    vecSpecies[ispec]->specMPI.patch_buff_index_send[idim][0].push_back( cuParticles.size()-1 );
                            //... without forgeting to add it to the exchange list for later cleaning...
	                    vecSpecies[ispec]->addPartInExchList(cuParticles.size()-1);
                            //... and removing it from receive buffer.
	                    (vecSpecies[ispec]->specMPI.patchVectorRecv[iDim][(iNeighbor+1)%2]).erase_particle(iPart);
	                    vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[iDim][(iNeighbor+1)%2]--;
                            check = 1;
	                }
                        //Other side of idim
	                else if ( (vecSpecies[ispec]->specMPI.patchVectorRecv[iDim][(iNeighbor+1)%2]).position(idim,iPart) >= max_local[idim]) { 
	            	    if (smpi->periods_[idim]==1 && Pcoordinates[idim] == params.number_of_patches[idim]-1) {
	                        (vecSpecies[ispec]->specMPI.patchVectorRecv[iDim][(iNeighbor+1)%2]).position(idim,iPart)     -= xmax[idim];
                            }
                            (vecSpecies[ispec]->specMPI.patchVectorRecv[iDim][(iNeighbor+1)%2]).cp_particle(iPart, cuParticles);
	                    vecSpecies[ispec]->specMPI.patch_buff_index_send[idim][1].push_back( cuParticles.size()-1 );
	                    vecSpecies[ispec]->addPartInExchList(cuParticles.size()-1);
                            (*cubmax)[(*cubmax).size()-1]++;
	                    (vecSpecies[ispec]->specMPI.patchVectorRecv[iDim][(iNeighbor+1)%2]).erase_particle(iPart);
	                    vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[iDim][(iNeighbor+1)%2]--;
                            check = 1;
	                }
                        idim++;
                    }
                }
                //We need to erase particles not in this domain from the received buffer
	        //for (int iPart=n_part_recv-1 ; iPart>=0; iPart-- ) {
	        //    if ( !(vecSpecies[ispec]->specMPI.patchVectorRecv[iDim][(iNeighbor+1)%2]).is_part_in_domain(iPart, this) ) {
	        //        (vecSpecies[ispec]->specMPI.patchVectorRecv[iDim][(iNeighbor+1)%2]).erase_particle(iPart);
	        //        vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[iDim][(iNeighbor+1)%2]--;
	        //    }
	        //}

            }//If not last dim for diagonal particles.

	} //If received something
    } //loop i Neighbor


    //        
    //    
    ///********************************************************************************/
    //// Delete Particles included in buff_send/buff_recv
    ///********************************************************************************/
    //cleanup_sent_particles(ispec, indexes_of_particles_to_exchange);
    //(*indexes_of_particles_to_exchange).clear();


    //// Delete useless Particles
    ////Theoretically, not even necessary to do anything as long you use bmax as the end of your iterator on particles.
    ////Nevertheless, you might want to free memory and have the actual number of particles
    ////really equal to the size of the vector. So we do:
    //cuParticles.erase_particle_trail((*cubmax).back());
    ////+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //La recopie finale doit se faire au traitement de la dernière dimension seulement !!
    if (iDim == ndim-1){
   
        //Evaluation of the necessary shift of all bins.
        for (unsigned int j=0; j<(*cubmax).size()+1 ;j++){
            shift[j]=0;
        }

        for (idim == 0; idim < ndim; idim++){
            if (idim==0) {
                shift[1] += vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[idim][0];//Particles coming from south all go to bin 0 and shift all the other bins.
                shift[(*cubmax).size()] += vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[idim][1];//Used only to count the total number of particles arrived.
            } else{
                for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
                    n_part_recv = vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[idim][iNeighbor];
                    for (unsigned int j=0; j<n_part_recv ;j++){
                    //We first evaluate how many particles arrive in each bin.
                	ii = int((vecSpecies[ispec]->specMPI.patchVectorRecv[idim][iNeighbor].position(0,j)-min_local[0])/dbin);//bin in which the particle goes.
                	shift[ii+1]++; // It makes the next bins shift.
                    }
                }
            }
        }


        //Must be done sequentially
        for (unsigned int j=1; j<(*cubmax).size()+1;j++){ //bin 0 is not shifted.Last element of shift stores total number of arriving particles.
            shift[j]+=shift[j-1];
        }
        //Make room for new particles
        cuParticles.initialize( cuParticles.size()+shift[(*cubmax).size()], cuParticles.dimension() );
            
        //Shift bins, must be done sequentially
        for (unsigned int j=(*cubmax).size()-1; j>=1; j--){
            n_particles = (*cubmax)[j]-(*cubmin)[j]; //Nbr of particle in this bin
            nmove = min(n_particles,shift[j]); //Nbr of particles to move
            lmove = max(n_particles,shift[j]); //How far particles must be shifted
            if (nmove>0) cuParticles.overwrite_part2D((*cubmin)[j], (*cubmin)[j]+lmove, nmove);
            (*cubmin)[j] += shift[j];
            (*cubmax)[j] += shift[j];
        }
            
        //Space has been made now to write the arriving particles into the correct bins
        //iDim == 0  is the easy case, when particles arrive either in first or last bin.
        for (idim == 0; idim < ndim; idim++){
            if (iDim==0) {
                for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
                    n_part_recv = vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[idim][iNeighbor];
                    if ( (neighbor_[0][iNeighbor]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
                	ii = iNeighbor*((*cubmax).size()-1);//0 if iNeighbor=0(particles coming from West) and (*cubmax).size()-1 otherwise.
                	vecSpecies[ispec]->specMPI.patchVectorRecv[idim][iNeighbor].overwrite_part2D(0, cuParticles,(*cubmax)[ii],n_part_recv);
                	(*cubmax)[ii] += n_part_recv ;
                    }
                }
            } else {
            //this is the difficult case, when particles can arrive in any bin.
                for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
                    n_part_recv = vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[idim][iNeighbor];
                    if ( (neighbor_[1][iNeighbor]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
                	for(unsigned int j=0; j<n_part_recv; j++){
                	    ii = int((vecSpecies[ispec]->specMPI.patchVectorRecv[iDim][iNeighbor].position(0,j)-min_local[0])/dbin);//bin in which the particle goes.
                	    vecSpecies[ispec]->specMPI.patchVectorRecv[idim][iNeighbor].overwrite_part2D(j, cuParticles,(*cubmax)[ii]);
                	    (*cubmax)[ii] ++ ;
                	}
                    }
                }
            }
        }
    }//End Recv_buffers ==> particles

    // Inject corner particles at the end of the list, update bmax
    //if (iDim==cuParticles.dimension()-1) cout << "Number of diag particles " << diagonalParticles.size() << endl;
    //for (int iPart = 0 ; iPart<diagonalParticles.size() ; iPart++) {
    //    diagonalParticles.cp_particle(iPart, cuParticles);
    //    (*indexes_of_particles_to_exchange).push_back(cuParticles.size()-1);
    //    (*cubmax)[(*cubmax).size()-1]++;
    //}


} // finalizeCommParticles(... iDim)


void Patch::createType( PicParams& params )
{
    int nx0 = params.n_space[0] + 1 + 2*params.oversize[0];
    int ny0 = params.n_space[1] + 1 + 2*params.oversize[1];
    unsigned int clrw = params.clrw;
    
    // MPI_Datatype ntype_[nDim][primDual][primDual]
    int nx, ny;
    int nline, ncol;

    int corner_nx, corner_ny;
    for (int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++) {
        nx = nx0 + ix_isPrim;
	corner_nx = 1 + 2*params.oversize[0] + ix_isPrim;
        for (int iy_isPrim=0 ; iy_isPrim<2 ; iy_isPrim++) {
            ny = ny0 + iy_isPrim;
	    corner_ny = 1 + 2*params.oversize[1] + iy_isPrim;

	    // Standard Type
            ntype_[0][ix_isPrim][iy_isPrim] = NULL;
            MPI_Type_contiguous(ny, MPI_DOUBLE, &(ntype_[0][ix_isPrim][iy_isPrim]));    //line
            MPI_Type_commit( &(ntype_[0][ix_isPrim][iy_isPrim]) );
            ntype_[1][ix_isPrim][iy_isPrim] = NULL;
            MPI_Type_vector(nx, 1, ny, MPI_DOUBLE, &(ntype_[1][ix_isPrim][iy_isPrim])); // column
            MPI_Type_commit( &(ntype_[1][ix_isPrim][iy_isPrim]) );
            ntype_[2][ix_isPrim][iy_isPrim] = NULL;
            MPI_Type_contiguous(ny*clrw, MPI_DOUBLE, &(ntype_[2][ix_isPrim][iy_isPrim]));   //clrw lines
            MPI_Type_commit( &(ntype_[2][ix_isPrim][iy_isPrim]) );

            ntypeSum_[0][ix_isPrim][iy_isPrim] = NULL;
            nline = 1 + 2*params.oversize[0] + ix_isPrim;
            MPI_Type_contiguous(nline, ntype_[0][ix_isPrim][iy_isPrim], &(ntypeSum_[0][ix_isPrim][iy_isPrim]));    //line
            MPI_Type_commit( &(ntypeSum_[0][ix_isPrim][iy_isPrim]) );
            ntypeSum_[1][ix_isPrim][iy_isPrim] = NULL;
            ncol  = 1 + 2*params.oversize[1] + iy_isPrim;
            MPI_Type_vector(nx, ncol, ny, MPI_DOUBLE, &(ntypeSum_[1][ix_isPrim][iy_isPrim])); // column
            MPI_Type_commit( &(ntypeSum_[1][ix_isPrim][iy_isPrim]) );

	    // Corner Types
            corner_ntype_[0][ix_isPrim][iy_isPrim] = NULL;
            MPI_Type_contiguous(corner_ny, MPI_DOUBLE, &(corner_ntype_[0][ix_isPrim][iy_isPrim]));    //line
	    MPI_Type_commit( &(corner_ntype_[0][ix_isPrim][iy_isPrim]) );
            corner_ntype_[1][ix_isPrim][iy_isPrim] = NULL;
            MPI_Type_vector(corner_nx, 1, ny, MPI_DOUBLE, &(corner_ntype_[1][ix_isPrim][iy_isPrim])); // column
	    MPI_Type_commit( &(corner_ntype_[1][ix_isPrim][iy_isPrim]) );
	    // clrw slide don't use corners
	    //corner_ntype_[2][ix_isPrim][iy_isPrim] = NULL;
            //MPI_Type_contiguous(ny*clrw, MPI_DOUBLE, &(ntype_[2][ix_isPrim][iy_isPrim]));   //clrw lines
	    //MPI_Type_commit( &(ntype_[2][ix_isPrim][iy_isPrim]) );

            corner_ntypeSum_[0][ix_isPrim][iy_isPrim] = NULL;
            nline = 1 + 2*params.oversize[0] + ix_isPrim;
	    MPI_Type_vector( nline, corner_ny, ny, MPI_DOUBLE, &(corner_ntypeSum_[0][ix_isPrim][iy_isPrim])); // column
	    MPI_Type_commit( &(corner_ntypeSum_[0][ix_isPrim][iy_isPrim]) );
            corner_ntypeSum_[1][ix_isPrim][iy_isPrim] = NULL;
	    ncol  = 1 + 2*params.oversize[1] + iy_isPrim;
	    MPI_Type_vector(corner_nx, ncol, ny, MPI_DOUBLE, &(corner_ntypeSum_[1][ix_isPrim][iy_isPrim])); // column
	    MPI_Type_commit( &(corner_ntypeSum_[1][ix_isPrim][iy_isPrim]) );

            
        }
    }
    
} //END createType

//Useless function ?
//void Patch::initSumRhoJ( ElectroMagn* EMfields, unsigned int diag_flag )
//{
//    // sum total charge density and currents
//
//    if (diag_flag) initSumField( EMfields->rho_ );
//    initSumField( EMfields->Jx_ );
//    initSumField( EMfields->Jy_ );
//    initSumField( EMfields->Jz_ );
//
//}

void Patch::initSumField( Field* field, int iDim )
{
    int patch_ndims_(2);
    int patch_nbNeighbors_(2);
    vector<unsigned int> patch_oversize(2,2);
    
    
    std::vector<unsigned int> n_elem = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field2D* f2D =  static_cast<Field2D*>(field);
   
    // Use a buffer per direction to exchange data before summing
    //Field2D buf[patch_ndims_][patch_nbNeighbors_];
    //Field2D corner_buf[patch_ndims_][patch_nbNeighbors_];
    // Size buffer is 2 oversize (1 inside & 1 outside of the current subdomain)
    std::vector<unsigned int> oversize2 = patch_oversize;
    oversize2[0] *= 2;
    oversize2[0] += 1 + f2D->isDual_[0];
    oversize2[1] *= 2;
    oversize2[1] += 1 + f2D->isDual_[1];
    
    //for (int iDim=0 ; iDim<patch_ndims_ ; iDim++) {
        for (int iNeighbor=0 ; iNeighbor<patch_nbNeighbors_ ; iNeighbor++) {
            std::vector<unsigned int> tmp(patch_ndims_,0);
            tmp[0] =    iDim  * n_elem[0] + (1-iDim) * oversize2[0];
            tmp[1] = (1-iDim) * n_elem[1] +    iDim  * oversize2[1];
            buf[iDim][iNeighbor].allocateDims( tmp );
        }
	//}
    /* for (int iDim=0 ; iDim<patch_ndims_ ; iDim++) {
        for (int iNeighbor=0 ; iNeighbor<patch_nbNeighbors_ ; iNeighbor++) {
            std::vector<unsigned int> tmp(patch_ndims_,0);
            tmp[0] = 1 + 2 * patch_oversize[0] + isDual[0];
            tmp[1] = 1 + 2 * patch_oversize[1] + isDual[1];
            corner_buf[iDim][iNeighbor].allocateDims( tmp );
        }
    }*/
     
    int istart, ix, iy;
    /********************************************************************************/
    // Send/Recv in a buffer data to sum
    /********************************************************************************/
    //for (int iDim=0 ; iDim<2 ; iDim++) {
        
	MPI_Datatype ntype = ntypeSum_[iDim][isDual[0]][isDual[1]];
	MPI_Request srequest[patch_ndims_][2];
	MPI_Request rrequest[patch_ndims_][2];
        
	for (int iNeighbor=0 ; iNeighbor<patch_nbNeighbors_ ; iNeighbor++) {
            
	    if ( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
		istart = iNeighbor * ( n_elem[iDim]- oversize2[iDim] ) + (1-iNeighbor) * ( 0 );
		ix = (1-iDim)*istart;
		iy =    iDim *istart;
		//int tag = buildtag( 3, hindex, neighbor_[iDim][iNeighbor]);
		int tag = buildtag( hindex, iDim, iNeighbor );
		//cout << hindex << " send to " << neighbor_[iDim][iNeighbor] << endl;
		//MPI_Isend( &(f2D->data_2D[ix][iy]), 1, ntype, 0, tag, MPI_COMM_SELF, &(f2D->specMPI.patch_srequest[iDim][iNeighbor]) );
		MPI_Isend( &(f2D->data_2D[ix][iy]), 1, ntype, MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &(f2D->specMPI.patch_srequest[iDim][iNeighbor]) );
	    } // END of Send
            
	    if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {
		int tmp_elem = (buf[iDim][(iNeighbor+1)%2]).dims_[0]*(buf[iDim][(iNeighbor+1)%2]).dims_[1];
		//int tag = buildtag( 3, neighbor_[iDim][(iNeighbor+1)%2], hindex);
		int tag = buildtag( neighbor_[iDim][(iNeighbor+1)%2], iDim, iNeighbor );
		//cout << hindex << " recv from " << neighbor_[iDim][(iNeighbor+1)%2] << " ; n_elements = " << tmp_elem << endl;
		//MPI_Irecv( &( (buf[iDim][(iNeighbor+1)%2]).data_2D[0][0] ), tmp_elem, MPI_DOUBLE, 0, tag, MPI_COMM_SELF, &(f2D->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]) );
		MPI_Irecv( &( (buf[iDim][(iNeighbor+1)%2]).data_2D[0][0] ), tmp_elem, MPI_DOUBLE, MPI_neighbor_[iDim][(iNeighbor+1)%2], tag, MPI_COMM_WORLD, &(f2D->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]) );
	    } // END of Recv
            
	} // END for iNeighbor
	//}

    /*for (int iDim=0 ; iDim<2 ; iDim++) {
	
	MPI_Datatype ntype = corner_ntypeSum_[0][isDual[0]][isDual[1]]; // 1st dimension useless
      
	for (int iNeighbor=0 ; iNeighbor<patch_nbNeighbors_ ; iNeighbor++) {
            
	    if ( is_a_corner_MPI_neighbor( iDim, iNeighbor ) ) {
		ix = iDim      * ( n_elem[0]     - oversize2[0]      ) + (1-iDim     ) * ( 0 );
		iy = iNeighbor * ( n_elem[1]- oversize2[1] ) + (1-iNeighbor) * ( 0 );
		int tag = buildtag( 3, hindex, corner_neighbor_[iDim][iNeighbor]);
		int tabsize(0);
		MPI_Type_size( ntype, &tabsize );
		//cout << hindex << " send in diagonal to " << corner_neighbor_[iDim][iNeighbor] << " from " << ix << " " << iy << " " << tabsize <<  endl;
		//MPI_Isend( &(f2D->data_2D[ix][iy]), 1, ntype, 0, tag, MPI_COMM_SELF, &(f2D->specMPI.corner_srequest[iDim][iNeighbor]) );
		MPI_Isend( &(f2D->data_2D[ix][iy]), 1, ntype, MPI_corner_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &(f2D->specMPI.corner_srequest[iDim][iNeighbor]) );
	    } // END of Send
            
	    if ( is_a_corner_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {
		int tmp_elem = (corner_buf[iDim][(iNeighbor+1)%2]).dims_[0]*(corner_buf[iDim][(iNeighbor+1)%2]).dims_[1];
		int tag = buildtag( 3, corner_neighbor_[iDim][(iNeighbor+1)%2], hindex);
		//cout << hindex << " recv from " << corner_neighbor_[iDim][(iNeighbor+1)%2] << " ; n_elements = " << tmp_elem << endl;
		//MPI_Irecv( &( (corner_buf[iDim][(iNeighbor+1)%2]).data_2D[0][0] ), tmp_elem, MPI_DOUBLE, 0, tag, MPI_COMM_SELF, &(f2D->specMPI.corner_rrequest[iDim][(iNeighbor+1)%2]) );
		MPI_Irecv( &( (corner_buf[iDim][(iNeighbor+1)%2]).data_2D[0][0] ), tmp_elem, MPI_DOUBLE, MPI_corner_neighbor_[iDim][(iNeighbor+1)%2], tag, MPI_COMM_WORLD, &(f2D->specMPI.corner_rrequest[iDim][(iNeighbor+1)%2]) );
	    } // END of Recv
            
	} // END for iNeighbor
    }*/

} // END initSumField

//Useless function ?
//void Patch::finalizeSumRhoJ( ElectroMagn* EMfields, unsigned int diag_flag )
//{
//    // sum total charge density and currents
//    if (diag_flag) finalizeSumField( EMfields->rho_ );
//    finalizeSumField( EMfields->Jx_ );
//    finalizeSumField( EMfields->Jy_ );
//    finalizeSumField( EMfields->Jz_ );
//
//}

void Patch::finalizeSumField( Field* field, int iDim )
{
    int patch_ndims_(2);
    int patch_nbNeighbors_(2);
    vector<unsigned int> patch_oversize(2,2);
    std::vector<unsigned int> n_elem = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field2D* f2D =  static_cast<Field2D*>(field);
   
    // Use a buffer per direction to exchange data before summing
    //Field2D buf[patch_ndims_][patch_nbNeighbors_];
    //Field2D corner_buf[patch_ndims_][patch_nbNeighbors_];
    // Size buffer is 2 oversize (1 inside & 1 outside of the current subdomain)
    std::vector<unsigned int> oversize2 = patch_oversize;
    oversize2[0] *= 2;
    oversize2[0] += 1 + f2D->isDual_[0];
    oversize2[1] *= 2;
    oversize2[1] += 1 + f2D->isDual_[1];
    
    int istart, ix, iy;
    /********************************************************************************/
    // Send/Recv in a buffer data to sum
    /********************************************************************************/

    //for (int iDim=0 ; iDim<2 ; iDim++) {
 
        MPI_Datatype ntype = ntypeSum_[iDim][isDual[0]][isDual[1]];
        MPI_Status sstat    [patch_ndims_][2];
        MPI_Status rstat    [patch_ndims_][2];
	
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
	    if ( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
		//cout << hindex << " is waiting for send at " << neighbor_[iDim][iNeighbor] << endl;
                MPI_Wait( &(f2D->specMPI.patch_srequest[iDim][iNeighbor]), &(sstat[iDim][iNeighbor]) );
            }
	    if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {
		//cout << hindex << " is waiting for recv from " << neighbor_[iDim][(iNeighbor+1)%2] << endl;	
                MPI_Wait( &(f2D->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]), &(rstat[iDim][(iNeighbor+1)%2]) );
            }
        }
	//}

    /*for (int iDim=0 ; iDim<2 ; iDim++) {
 
        MPI_Datatype ntype = corner_ntypeSum_[0][isDual[0]][isDual[1]];; // 1st dimension useless
        MPI_Status corner_sstat    [patch_ndims_][2];
        MPI_Status corner_rstat    [patch_ndims_][2];
	
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
	    if ( is_a_corner_MPI_neighbor( iDim, iNeighbor ) ) {
		//cout << hindex << " is waiting for corner send at " << corner_neighbor_[iDim][iNeighbor] << endl;
                MPI_Wait( &(f2D->specMPI.corner_srequest[iDim][iNeighbor]), &(corner_sstat[iDim][iNeighbor]) );
		//cout << hindex << " is waiting for corner send at " << corner_neighbor_[iDim][iNeighbor] << " ACHIEVED" << endl;
            }
	    if ( is_a_corner_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {
		//cout << hindex << " is waiting for corner recv from " << corner_neighbor_[iDim][(iNeighbor+1)%2] << endl;	
                MPI_Wait( &(f2D->specMPI.corner_rrequest[iDim][(iNeighbor+1)%2]), &(corner_rstat[iDim][(iNeighbor+1)%2]) );
            }
        }
    }*/

    /********************************************************************************/
    // Sum data on each process, same operation on both side
    /********************************************************************************/
    //for (int iDim=0 ; iDim<2 ; iDim++) {
       
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
            istart = ( (iNeighbor+1)%2 ) * ( n_elem[iDim]- oversize2[iDim] ) + (1-(iNeighbor+1)%2) * ( 0 );
            int ix0 = (1-iDim)*istart;
            int iy0 =    iDim *istart;
	    if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {
                for (unsigned int ix=0 ; ix< (buf[iDim][(iNeighbor+1)%2]).dims_[0] ; ix++) {
                    for (unsigned int iy=0 ; iy< (buf[iDim][(iNeighbor+1)%2]).dims_[1] ; iy++)
                        f2D->data_2D[ix0+ix][iy0+iy] += (buf[iDim][(iNeighbor+1)%2])(ix,iy);
                }
            } // END if
            
        } // END for iNeighbor
	//}
        
    /*for (int iDim=0 ; iDim<2 ; iDim++) {
       
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
	    int ix0 = iDim      * ( n_elem[0]     - oversize2[0]      ) + (1-iDim     ) * ( 0 );
	    int iy0 = iNeighbor * ( n_elem[1]- oversize2[1] ) + (1-iNeighbor) * ( 0 );
	    if ( is_a_corner_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {
                for (unsigned int ix=0 ; ix< (corner_buf[iDim][(iNeighbor+1)%2]).dims_[0] ; ix++) {
                    for (unsigned int iy=0 ; iy< (corner_buf[iDim][(iNeighbor+1)%2]).dims_[1] ; iy++)
                        f2D->data_2D[ix0+ix][iy0+iy] += (corner_buf[iDim][(iNeighbor+1)%2])(ix,iy);
                }
            } // END if
            
        } // END for iNeighbor
    } */      

    //for (int iDim=0 ; iDim<patch_ndims_ ; iDim++) {
        for (int iNeighbor=0 ; iNeighbor<patch_nbNeighbors_ ; iNeighbor++) {
	    buf[iDim][iNeighbor].deallocateDims();
        }
	//}
    /*for (int iDim=0 ; iDim<patch_ndims_ ; iDim++) {
        for (int iNeighbor=0 ; iNeighbor<patch_nbNeighbors_ ; iNeighbor++) {
	    corner_buf[iDim][iNeighbor].deallocateDims();
        }
    }*/

} // END finalizeSumField

void Patch::initExchange( Field* field )
{
    int patch_ndims_(2);
    int patch_nbNeighbors_(2);
    vector<unsigned int> patch_oversize(2,2);

    std::vector<unsigned int> n_elem   = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field2D* f2D =  static_cast<Field2D*>(field);

    int istart, ix, iy;

    // Loop over dimField
    for (int iDim=0 ; iDim<patch_ndims_ ; iDim++) {

        MPI_Datatype ntype = ntype_[iDim][isDual[0]][isDual[1]];
        for (int iNeighbor=0 ; iNeighbor<patch_nbNeighbors_ ; iNeighbor++) {

	    if ( is_a_MPI_neighbor( iDim, iNeighbor ) ) {

                istart = iNeighbor * ( n_elem[iDim]- (2*patch_oversize[iDim]+1+isDual[iDim]) ) + (1-iNeighbor) * ( 2*patch_oversize[iDim] + isDual[iDim] );
                ix = (1-iDim)*istart;
                iy =    iDim *istart;
		//int tag = buildtag( 4, hindex, neighbor_[iDim][iNeighbor]);
		int tag = buildtag( hindex, iDim, iNeighbor );
                //MPI_Isend( &(f2D->data_2D[ix][iy]), 1, ntype, 0, tag, MPI_COMM_SELF, &(f2D->specMPI.patch_srequest[iDim][iNeighbor]) );
                MPI_Isend( &(f2D->data_2D[ix][iy]), 1, ntype, MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &(f2D->specMPI.patch_srequest[iDim][iNeighbor]) );

            } // END of Send

	    if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {

                istart = ( (iNeighbor+1)%2 ) * ( n_elem[iDim] - 1 ) + (1-(iNeighbor+1)%2) * ( 0 )  ;
                ix = (1-iDim)*istart;
                iy =    iDim *istart;
 		//int tag = buildtag( 4, neighbor_[iDim][(iNeighbor+1)%2], hindex);
 		int tag = buildtag( neighbor_[iDim][(iNeighbor+1)%2], iDim, iNeighbor );
		//MPI_Irecv( &(f2D->data_2D[ix][iy]), 1, ntype, 0, tag, MPI_COMM_SELF, &(f2D->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]));
		MPI_Irecv( &(f2D->data_2D[ix][iy]), 1, ntype, MPI_neighbor_[iDim][(iNeighbor+1)%2], tag, MPI_COMM_WORLD, &(f2D->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]));

            } // END of Recv

        } // END for iNeighbor

    } // END for iDim

}
void Patch::initExchange( Field* field, int iDim )
{
    int patch_ndims_(2);
    int patch_nbNeighbors_(2);
    vector<unsigned int> patch_oversize(2,2);

    std::vector<unsigned int> n_elem   = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field2D* f2D =  static_cast<Field2D*>(field);

    int istart, ix, iy;

    MPI_Datatype ntype = ntype_[iDim][isDual[0]][isDual[1]];
    for (int iNeighbor=0 ; iNeighbor<patch_nbNeighbors_ ; iNeighbor++) {

	if ( is_a_MPI_neighbor( iDim, iNeighbor ) ) {

	    istart = iNeighbor * ( n_elem[iDim]- (2*patch_oversize[iDim]+1+isDual[iDim]) ) + (1-iNeighbor) * ( 2*patch_oversize[iDim] + isDual[iDim] );
	    ix = (1-iDim)*istart;
	    iy =    iDim *istart;
	    int tag = buildtag( hindex, iDim, iNeighbor );
	    //MPI_Isend( &(f2D->data_2D[ix][iy]), 1, ntype, 0, tag, MPI_COMM_SELF, &(f2D->specMPI.patch_srequest[iDim][iNeighbor]) );
	    MPI_Isend( &(f2D->data_2D[ix][iy]), 1, ntype, MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &(f2D->specMPI.patch_srequest[iDim][iNeighbor]) );

	} // END of Send

	if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {

	    istart = ( (iNeighbor+1)%2 ) * ( n_elem[iDim] - 1 ) + (1-(iNeighbor+1)%2) * ( 0 )  ;
	    ix = (1-iDim)*istart;
	    iy =    iDim *istart;
	    int tag = buildtag( neighbor_[iDim][(iNeighbor+1)%2], iDim, iNeighbor );
	    //MPI_Irecv( &(f2D->data_2D[ix][iy]), 1, ntype, 0, tag, MPI_COMM_SELF, &(f2D->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]));
	    MPI_Irecv( &(f2D->data_2D[ix][iy]), 1, ntype, MPI_neighbor_[iDim][(iNeighbor+1)%2], tag, MPI_COMM_WORLD, &(f2D->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]));

	} // END of Recv

    } // END for iNeighbor


}


void Patch::finalizeExchange( Field* field )
{
    int patch_ndims_(2);
    int patch_nbNeighbors_(2);

    Field2D* f2D =  static_cast<Field2D*>(field);

    MPI_Status sstat    [patch_ndims_][2];
    MPI_Status rstat    [patch_ndims_][2];

    // Loop over dimField
    for (int iDim=0 ; iDim<patch_ndims_ ; iDim++) {

        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
	    if ( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
                MPI_Wait( &(f2D->specMPI.patch_srequest[iDim][iNeighbor]), &(sstat[iDim][iNeighbor]) );
            }
 	    if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {
               MPI_Wait( &(f2D->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]), &(rstat[iDim][(iNeighbor+1)%2]) );
            }
        }

    } // END for iDim
}

void Patch::finalizeExchange( Field* field, int iDim )
{
    int patch_ndims_(2);
    int patch_nbNeighbors_(2);

    Field2D* f2D =  static_cast<Field2D*>(field);

    MPI_Status sstat    [patch_ndims_][2];
    MPI_Status rstat    [patch_ndims_][2];

    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
	if ( is_a_MPI_neighbor( iDim, iNeighbor ) ) {
	    MPI_Wait( &(f2D->specMPI.patch_srequest[iDim][iNeighbor]), &(sstat[iDim][iNeighbor]) );
	}
	if ( is_a_MPI_neighbor( iDim, (iNeighbor+1)%2 ) ) {
	    MPI_Wait( &(f2D->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]), &(rstat[iDim][(iNeighbor+1)%2]) );
	}
    }

}

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
		cuParticles.overwrite_part2D((*cubmax)[ibin]-1, iPart );
		(*cubmax)[ibin]--;
		ii--;
		iPart = (*indexes_of_particles_to_exchange)[ii];
	    }
	    if (iPart >= (*cubmin)[ibin] && iPart < (*cubmax)[ibin]) { //On traite la dernière particule (qui peut aussi etre la premiere)
		cuParticles.overwrite_part2D((*cubmax)[ibin]-1, iPart );
		(*cubmax)[ibin]--;
	    }
	}
    }


    //Shift the bins in memory
    //Warning: this loop must be executed sequentially. Do not use openMP here.
    for (int unsigned ibin = 1 ; ibin < (*cubmax).size() ; ibin++ ) { //First bin don't need to be shifted
	ii = (*cubmin)[ibin]-(*cubmax)[ibin-1]; // Shift the bin in memory by ii slots.
	iPart = min(ii,(*cubmax)[ibin]-(*cubmin)[ibin]); // Number of particles we have to shift = min (Nshift, Nparticle in the bin)
	if(iPart > 0) cuParticles.overwrite_part2D((*cubmax)[ibin]-iPart,(*cubmax)[ibin-1],iPart);
	(*cubmax)[ibin] -= ii;
	(*cubmin)[ibin] = (*cubmax)[ibin-1];
    }

}
