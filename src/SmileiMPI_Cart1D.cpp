#include "SmileiMPI_Cart1D.h"

#include "Species.h"
#include "ParticleFactory.h"

#include "ElectroMagn.h"
#include "Field1D.h"

#include "Tools.h" 

#include <string>
#include <cmath>

#include "Field2D.h"
#include <string.h>

using namespace std;

SmileiMPI_Cart1D::SmileiMPI_Cart1D( int* argc, char*** argv )
	: SmileiMPI( argc, argv )
{
}

SmileiMPI_Cart1D::SmileiMPI_Cart1D( SmileiMPI* smpi)
	: SmileiMPI( smpi )
{
	ndims_ = 1;
	number_of_procs  = new int(ndims_);
	coords_  = new int(ndims_);
	periods_  = new int(ndims_);
	reorder_ = 0;

	nbNeighbors_ = 2;

	for (int i=0 ; i<ndims_ ; i++) periods_[i] = 0;
	for (int i=0 ; i<ndims_ ; i++) coords_[i] = 0;
	for (int i=0 ; i<ndims_ ; i++) number_of_procs[i] = 1;

	for (int iDim=0 ; iDim<ndims_ ; iDim++) {
		for (int i=0 ; i<nbNeighbors_ ; i++) {
			neighbor_[iDim][i] = MPI_PROC_NULL;
			buff_index_send[iDim][i].resize(0);
			buff_index_recv_sz[iDim][i] = 0;
		}
	}

}

SmileiMPI_Cart1D::~SmileiMPI_Cart1D()
{
	delete number_of_procs;
	delete periods_;
	delete coords_;

	if ( SMILEI_COMM_1D != MPI_COMM_NULL) MPI_Comm_free(&SMILEI_COMM_1D);

}

void SmileiMPI_Cart1D::createTopology(PicParams& params)
{
    for (unsigned int i=0 ; i<params.nDim_field ; i++)
		params.n_space_global[i] = round(params.res_space[i]*params.sim_length[i]/(2.0*M_PI));

    number_of_procs[0] = smilei_sz;
    
	MPI_Cart_create( SMILEI_COMM_WORLD, ndims_, number_of_procs, periods_, reorder_, &SMILEI_COMM_1D );
	MPI_Cart_coords( SMILEI_COMM_1D, smilei_rk, ndims_, coords_ );
    
    
	for (int iDim=0 ; iDim<ndims_ ; iDim++) {
		MPI_Cart_shift( SMILEI_COMM_1D, iDim, 1, &(neighbor_[iDim][0]), &(neighbor_[iDim][1]) );
		//PMESSAGE ( 0, smilei_rk, "Neighbors of process in direction " << iDim << " : " << neighbor_[iDim][0] << " - " << neighbor_[iDim][1]  );
	}
    
    
	for (unsigned int i=0 ; i<params.nDim_field ; i++) {
        
		params.n_space[i] = params.n_space_global[i] / number_of_procs[i];
        
		n_space_global[i] = params.n_space_global[i];
		oversize[i] = params.oversize[i] = 2;
		cell_starting_global_index[i] = coords_[i]*(params.n_space_global[i] / number_of_procs[i]);
        
        
		if ( number_of_procs[i]*params.n_space[i] != params.n_space_global[i] ) {
			// Correction on the last MPI process of the direction to use the wished number of cells
			if (coords_[i]==number_of_procs[i]-1) {
				params.n_space[i] = params.n_space_global[i] - params.n_space[i]*(number_of_procs[i]-1);
			}
		}
		// min/max_local : describe local domain in which particles cat be moved
		//                 different from domain on which E, B, J are defined
		min_local[i] = (cell_starting_global_index[i]                  )*params.cell_length[i];
		max_local[i] = (cell_starting_global_index[i]+params.n_space[i])*params.cell_length[i];
		//PMESSAGE( 0, smilei_rk, "min_local / mac_local on " << smilei_rk << " = " << min_local[i] << " / " << max_local[i] << " selon la direction " << i );
        
		cell_starting_global_index[i] -= params.oversize[i];
        
	}
    
    
    
	MESSAGE( "n_space / rank " << smilei_rk << " = " << params.n_space[0]  );

	extrem_ranks[0][0] = 0;
	int rank_min =  0;
   	if (coords_[0] == 0) {
		rank_min = smilei_rk;
	} 
	MPI_Allreduce(&rank_min, &extrem_ranks[0][0], 1, MPI_INT, MPI_SUM, SMILEI_COMM_1D);
	extrem_ranks[0][1] = 0;
	int rank_max = 0;
	if (coords_[0]==number_of_procs[0]-1) {
		rank_max = smilei_rk;
	}
	MPI_Allreduce(&rank_max, &extrem_ranks[0][1], 1, MPI_INT, MPI_SUM, SMILEI_COMM_1D);

	//cout << extrem_ranks[0][0] << " " << extrem_ranks[0][1] << endl;
    
    
/*	number_of_procs[0] = smilei_sz;

	MPI_Cart_create( SMILEI_COMM_WORLD, ndims_, number_of_procs, periods_, reorder_, &SMILEI_COMM_1D );
	MPI_Cart_coords( SMILEI_COMM_1D, smilei_rk, ndims_, coords_ );

	// neighbor_[0][0]  |  Current process  |  neighbor_[0][1] //
	MPI_Cart_shift( SMILEI_COMM_1D, 0, 1, &(neighbor_[0][0]), &(neighbor_[0][1]) );
	PMESSAGE ( 0, smilei_rk, "Neighbors of process : " << neighbor_[0][0] << " - " << neighbor_[0][1]  );


	for (unsigned int i=0 ; i<params.nDim_field ; i++) {

		params.n_space[i] = params.n_space_global[i] / number_of_procs[i];
//		if ( number_of_procs[i]*params.n_space[i] != params.n_space_global[i] ) {
//			//WARNING( "Domain splitting does not match to the global domain" );
//			if (coords_[i]==number_of_procs[i]-1) {
//				params.n_space[i] = params.n_space_global[i] - params.n_space[i]*(number_of_procs[i]-1);
//			}
//		}

		n_space_global[i] = params.n_space_global[i];
		oversize[i] = params.oversize[i] = 2;
		//! \todo{replace cell_starting_global_index compute by a most sophisticated or input data}
		cell_starting_global_index[i] = coords_[i]*params.n_space[i];
		// min/max_local : describe local domain in which particles cat be moved
		//                 different from domain on which E, B, J are defined
		min_local[i] = (cell_starting_global_index[i]                  )*params.cell_length[i];
		max_local[i] = (cell_starting_global_index[i]+params.n_space[i])*params.cell_length[i];
		cell_starting_global_index[i] -= params.oversize[i];

		if ( number_of_procs[i]*params.n_space[i] != params.n_space_global[i] ) {
			//WARNING( "Domain splitting does not match to the global domain" );
			if (coords_[i]==number_of_procs[i]-1) {
				params.n_space[i] = params.n_space_global[i] - params.n_space[i]*(number_of_procs[i]-1);
			}
		}

	}*/
}

void SmileiMPI_Cart1D::exchangeParticles(Species* species, int ispec, PicParams* params)
{
	std::vector<Particle*>* cuParticles = &species->particles;
        std::vector<int>* cubmin = &species->bmin;
        std::vector<int>* cubmax = &species->bmax;

        MPI_Status Stat;
        int n_particles;
        double dbin;

	/********************************************************************************/
	// Build lists of indexes of particle to exchange per neighbor
	// Computed from indexes_of_particles_to_exchange computed during particles' BC
	/********************************************************************************/
	int n_part_send = indexes_of_particles_to_exchange.size();
	int n_part_recv;
        
	int ii, iPart;
	for (int i=0 ; i<n_part_send ; i++) {
		iPart = indexes_of_particles_to_exchange[i];
		if      ( (*cuParticles)[iPart]->position(0) < min_local[0]){
			buff_index_send[0][0].push_back( indexes_of_particles_to_exchange[i] );
                }
		else if ( (*cuParticles)[iPart]->position(0) >= max_local[0]){
			buff_index_send[0][1].push_back( indexes_of_particles_to_exchange[i] );
                }
	} // END for iPart = f(i)


	/********************************************************************************/
	// Exchange particles
	/********************************************************************************/
	// Loop over neighbors in a direction

	MPI_Status stat[2];
	MPI_Request request[2];
	// iDim = 0
	// Send to neighbor_[0][iNeighbor] / Recv from neighbor_[0][(iNeighbor+1)%2] :
	// MPI_COMM_SIZE = 2 :  neighbor_[0][0]  |  Current process  |  neighbor_[0][1]
	// Rank = 0 : iNeighbor = 0 : neighbor_[0][0] = NONE : neighbor_[0][(0+1)%2 = 1
	//            iNeighbor = 1 : neighbor_[0][1] = 1    : neighbor_[0][(1+1)%2 = NONE
	// Rank = 1 : iNeighbor = 0 : neighbor_[0][0] = 0    : neighbor_[0][(0+1)%2 = NONE
	//            iNeighbor = 1 : neighbor_[0][1] = NONE : neighbor_[0][(1+1)%2 = 0

	///********************************************************************************/
	//// Exchange number of particles to exchange to establish or not a communication
	///********************************************************************************/

        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
	    n_part_send = buff_index_send[0][iNeighbor].size();
            if ( (neighbor_[0][0]!=MPI_PROC_NULL) && (neighbor_[0][1]!=MPI_PROC_NULL) ){
               //Send-receive
                 MPI_Sendrecv( &n_part_send, 1, MPI_INT, neighbor_[0][iNeighbor], 0, &buff_index_recv_sz[0][(iNeighbor+1)%2], 1, MPI_INT, neighbor_[0][(iNeighbor+1)%2], 0, SMILEI_COMM_1D,&Stat);
            } else if (neighbor_[0][iNeighbor]!=MPI_PROC_NULL){
                 //Send
                 MPI_Send( &n_part_send, 1, MPI_INT, neighbor_[0][iNeighbor], 0, SMILEI_COMM_1D);
            } else if (neighbor_[0][(iNeighbor+1)%2]!=MPI_PROC_NULL){
                 //Receive
                 MPI_Recv( &buff_index_recv_sz[0][(iNeighbor+1)%2], 1, MPI_INT, neighbor_[0][(iNeighbor+1)%2], 0, SMILEI_COMM_1D, &Stat);
            }
        }     

	/********************************************************************************/
	// Define buffers to exchange buff_index_send[iDim][iNeighbor].size();
	/********************************************************************************/	
	//! \todo Define this as a main parameter for the code so that it needs not be defined all the time
	int part_mem_size=(2*params->nDim_particle+3+1)*sizeof(double)+sizeof(short);
	
	std::vector<unsigned int> partSize(2,0);
	partSize[0] = 0; // Number of particle exchanged per direction
	//! \todo{7 replaced by number of properties by particle, short (charge) managed as a double }
	partSize[1] = 7;
	Field2D partArrayRecv[2];
	Field2D partArraySend[2];

	/********************************************************************************/
	// Proceed to effective Particles' communications
	/********************************************************************************/

        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
	    n_part_send = buff_index_send[0][iNeighbor].size();
	    n_part_recv = buff_index_recv_sz[0][(iNeighbor+1)%2];
            if ( (neighbor_[0][0]!=MPI_PROC_NULL) && (neighbor_[0][1]!=MPI_PROC_NULL) && (n_part_send!=0) && (n_part_recv!=0) ){
               //Send-receive
        	 partSize[0] = n_part_send;
		 partArraySend[iNeighbor].allocateDims(partSize);
		 for (int iPart=0 ; iPart<n_part_send ; iPart++) {
		     memcpy(&(partArraySend[iNeighbor](iPart,0)), &((*cuParticles)[ buff_index_send[0][iNeighbor][iPart] ]->position(0)), 7*sizeof(double) );
                 }
                 partSize[0] = n_part_recv;
		 partArrayRecv[(iNeighbor+1)%2].allocateDims(partSize);
	         MPI_Sendrecv( &(partArraySend[iNeighbor](0,0)), (int)7*n_part_send, MPI_DOUBLE, neighbor_[0][iNeighbor], 0, &(partArrayRecv[(iNeighbor+1)%2](0,0)), (int)7*n_part_recv, MPI_DOUBLE,  neighbor_[0][(iNeighbor+1)%2], 0, SMILEI_COMM_1D, &Stat);
            } else if ( (neighbor_[0][iNeighbor]!=MPI_PROC_NULL) && (n_part_send!=0) ){
                 //Send
        	 partSize[0] = n_part_send;
		 partArraySend[iNeighbor].allocateDims(partSize);
		 for (int iPart=0 ; iPart<n_part_send ; iPart++) {
		     memcpy(&(partArraySend[iNeighbor](iPart,0)), &((*cuParticles)[ buff_index_send[0][iNeighbor][iPart] ]->position(0)), 7*sizeof(double) );
                 }
	         MPI_Send( &(partArraySend[iNeighbor](0,0)), (int)7*n_part_send, MPI_DOUBLE, neighbor_[0][iNeighbor], 0, SMILEI_COMM_1D);
            } else if ( (neighbor_[0][(iNeighbor+1)%2]!=MPI_PROC_NULL) && (n_part_recv!=0) ){
                 //Receive
                 partSize[0] = n_part_recv;
		 partArrayRecv[(iNeighbor+1)%2].allocateDims(partSize);
		 MPI_Recv( &(partArrayRecv[(iNeighbor+1)%2](0,0)), (int)7*n_part_recv, MPI_DOUBLE,  neighbor_[0][(iNeighbor+1)%2], 0, SMILEI_COMM_1D, &Stat );
            }
        }  

       /********************************************************************************/
        // Delete Particles included in buff_send/buff_recv
        /********************************************************************************/
        // Push lost particles at the end of bins
        //! \todo For loop on bins, can use openMP here.
        for (unsigned int ibin = 0 ; ibin < (*cubmax).size() ; ibin++ ) {
            //cout << "bounds " << (*cubmin)[ibin] << " " << (*cubmax)[ibin] << endl;
            ii = indexes_of_particles_to_exchange.size()-1;
            if (ii >= 0){ // Push lost particles to the end of the bin
                iPart = indexes_of_particles_to_exchange[ii];
                while (iPart >= (*cubmax)[ibin] && ii > 0) {
                    ii--;
                    iPart = indexes_of_particles_to_exchange[ii];
                }
                while (iPart == (*cubmax)[ibin]-1 && ii > 0) {
                    (*cubmax)[ibin]--;
                    ii--;
                    iPart = indexes_of_particles_to_exchange[ii];
                }
                while (iPart >= (*cubmin)[ibin] && ii > 0) {
                    //!\todo swap can be switched to a simple overwrite with a single memcpy
                    species->swap_part( (*cuParticles)[iPart] , (*cuParticles)[(*cubmax)[ibin]-1] );
                    (*cubmax)[ibin]--;
                    ii--;
                    iPart = indexes_of_particles_to_exchange[ii];
                }
                if (iPart >= (*cubmin)[ibin] && iPart < (*cubmax)[ibin]) { //On traite la dernière particule (qui peut aussi etre la premiere)
                    //!\todo swap can be switched to a simple overwrite with a single memcpy
                    species->swap_part( (*cuParticles)[iPart] , (*cuParticles)[(*cubmax)[ibin]-1] );
                    (*cubmax)[ibin]--;
                }
                            }
        }
        //Shift the bins in memory
        //Warning: this loop must be executed sequentially. Do not use openMP here.
        //Here we can assume all particles are contiguous in memory and shift a full block of particle in only 1 memcpy.
        //If they are not, it must be decomposed into iPart memcpy.
        for (int unsigned ibin = 1 ; ibin < (*cubmax).size() ; ibin++ ) { //First bin don't need to be shifted
            ii = (*cubmin)[ibin]-(*cubmax)[ibin-1]; // Shift the bin in memory by ii slots. 
            iPart = min(ii,(*cubmax)[ibin]-(*cubmin)[ibin]); // Number of particles we have to shift = min (Nshift, Nparticle in the bin)
            //coniguous version:----------------------------------------------------------------------------------------------------------------------------
            //memcpy( &( ((*cuParticles)[(*cubmax)[ibin-1]])->position(0) ), &( ((*cuParticles)[(*cubmax)[ibin]-iPart])->position(0) ) , iPart*part_mem_size );
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            //non contiguous version:----------------------------------------------------------------------------------------------------------------------
            for (int i =0; i < iPart ; i++){
                memcpy( &( ((*cuParticles)[(*cubmax)[ibin-1]+i])->position(0) ), &( ((*cuParticles)[(*cubmax)[ibin]-i-1])->position(0) ) , part_mem_size );
            }
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            (*cubmax)[ibin] -= ii;
            (*cubmin)[ibin] = (*cubmax)[ibin-1];
        }

	// Delete useless Particles
        //contiguous version:------------------------------------------------------------
        //Not even necessary to do anything as long you use bmax as the end of your iterator on particles.
        //Nevertheless, you might want to free memory with something like:
        //cuParticles->erase( (*cuParticles)[(*cubmax).back()], cuParticles->end()-1 );
        //You need to erase all particles from cubmax.back() and beyond.
        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //non contiguous version:erase one by one ----------------------------------------
	n_particles = species->getNbrOfParticles();
	for (int i= n_particles-1; i >= (*cubmax).back() ; i--) {
		(*cuParticles)[i]->~Particle();
		//cuParticles->erase( cuParticles->begin() + i );
		(*cuParticles).pop_back();
	}
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


        /********************************************************************************/
        // Clean lists of indexes of particle to exchange per neighbor
        /********************************************************************************/
        for (int i=0 ; i<nbNeighbors_ ; i++)
                buff_index_send[0][i].clear();
        /********************************************************************************/
        // Copy newly arrived particles back to the vector
        /********************************************************************************/
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {

                n_part_recv = buff_index_recv_sz[0][(iNeighbor+1)%2];
                if ( (neighbor_[0][(iNeighbor+1)%2]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
                    if (iNeighbor == 0){ // Copy particles coming from the right at the end of Particles Array
                        n_particles = species->getNbrOfParticles();
                        cuParticles->resize( n_particles + n_part_recv );
                        for (int iPart=0 ; iPart<n_part_recv; iPart++ ) {
                                (*cuParticles)[n_particles+iPart] = ParticleFactory::create(params, ispec);
                                memcpy( &( ((*cuParticles)[n_particles+iPart])->position(0) ), &(partArrayRecv[(iNeighbor+1)%2](iPart,0)), part_mem_size );
                                //cout << ((*cuParticles)[n_particles+iPart])->position(0) << endl;
                        }
                        (*cubmax)[(*cubmax).size()-1] += n_part_recv ;
                    } else {// Copy particles coming from the left at the beginning of Particles Array
                        // This is extremely ugly and WILL have to be changed once all particles are contiguous in memory.
                        for (int iPart=0 ; iPart<n_part_recv; iPart++ ) {
                            //cout << "copy iPart "<< iPart <<" from the left" << endl;
                            (*cuParticles).insert((*cuParticles).begin(), ParticleFactory::create(params, ispec) );
                            memcpy( &( ((*cuParticles)[0])->position(0) ), &(partArrayRecv[(iNeighbor+1)%2](iPart,0)), part_mem_size );
                                //cout << ((*cuParticles)[0])->position(0) << endl;
                        }
                        (*cubmax)[0] += n_part_recv ;
                        for (unsigned int ibin=1 ; ibin < (*cubmax).size() ; ibin++ ) {
                            (*cubmax)[ibin] += n_part_recv ;
                            (*cubmin)[ibin] = (*cubmax)[ibin-1] ;
                        }
                    }
                }
        }

#ifdef _OLD
	/********************************************************************************/
	// Copy newly arrived particles back to the vector
	/********************************************************************************/
	for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {

		n_part_recv = buff_index_recv_sz[0][(iNeighbor+1)%2];
		
		if ( (neighbor_[0][(iNeighbor+1)%2]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
			int n_particles = species->getNbrOfParticles();
		
			cuParticles->resize( n_particles + n_part_recv );
			for (int iPart=0 ; iPart<n_part_recv; iPart++ ) {
				(*cuParticles)[n_particles+iPart] = ParticleFactory::create(params, ispec);
				memcpy( &( ((*cuParticles)[n_particles+iPart])->position(0) ), &(partArrayRecv[(iNeighbor+1)%2](iPart,0)), 7*sizeof(double) );
			}
		}
	}
	/********************************************************************************/
	// Clean lists of indexes of particle to exchange per neighbor
	/********************************************************************************/
	for (int i=0 ; i<nbNeighbors_ ; i++)
		buff_index_send[0][i].clear();

	/********************************************************************************/
	// Delete Particles included in buff_send/buff_recv
	/********************************************************************************/
	n_part_send = indexes_of_particles_to_exchange.size();
	for (int i=n_part_send-1 ; i>=0 ; i--) {
		iPart = indexes_of_particles_to_exchange[i];
		(*cuParticles)[iPart]->~Particle();
		cuParticles->erase( cuParticles->begin() + iPart );
	} // END for iPart = f(i)
#endif
	//DEBUG( 2, "\tProcess " << smilei_rk << " : " << species->getNbrOfParticles() << " Particles of species " << ispec );
} // END exchangeParticles

void SmileiMPI_Cart1D::sumField( Field* field )
{
	std::vector<unsigned int> n_elem = field->dims_;
	Field1D* f1D =  static_cast<Field1D*>(field);

	// Use a buffer per direction to exchange data before summing
	Field1D buf[ nbNeighbors_ ];
	// Size buffer is 2 oversize (1 inside & 1 outside of the current subdomain)
	std::vector<unsigned int> oversize2 = oversize;
	oversize2[0] *= 2; oversize2[0] += 1 + f1D->isPrimal_[0];
	for (int i=0;i<nbNeighbors_ ;i++)  buf[i].allocateDims( oversize2 );

	// istart store in the first part starting index of data to send, then the starting index of data to write in
	// Send point of vue : istart =           iNeighbor * ( n_elem[0]- 2*oversize[0] ) + (1-iNeighbor)       * ( 0 );
	// Rank = 0 : iNeighbor = 0 : send - neighbor_[0][0] = NONE
	//            iNeighbor = 1 : send - neighbor_[0][1] = 1 / istart = ( n_elem[0]- 2*oversize[0] )
	// Rank = 1 : iNeighbor = 0 : send - neighbor_[0][0] = 0 / istart = 0
	//            iNeighbor = 1 : send - neighbor_[0][1] = NONE
	// Recv point of vue : istart = ( (iNeighbor+1)%2 ) * ( n_elem[0]- 2*oversize[0] ) + (1-(iNeighbor+1)%2) * ( 0 );
	// Rank = 0 : iNeighbor = 0 : recv - neighbor_[0][1] = 1 / istart = ( n_elem[0]- 2*oversize[0] )
	//            iNeighbor = 1 : recv - neighbor_[0][0] = NONE
	// Rank = 1 : iNeighbor = 0 : recv - neighbor_[0][1] = NONE
	//            iNeighbor = 1 : recv - neighbor_[0][0] = 0 / istart = 0
	int istart;

	MPI_Status sstat[2];
	MPI_Status rstat[2];
	MPI_Request srequest[2];
	MPI_Request rrequest[2];
	/********************************************************************************/
	// Send/Recv in a buffer data to sum
	/********************************************************************************/
	// Loop over neighbors in a direction
	// Send to neighbor_[0][iNeighbor] / Recv from neighbor_[0][(iNeighbor+1)%2] :
	// See in exchangeParticles()
	for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {

		if (neighbor_[0][iNeighbor]!=MPI_PROC_NULL) {
			istart = iNeighbor * ( n_elem[0]- oversize2[0] ) + (1-iNeighbor) * ( 0 );
			MPI_Isend( &(f1D->data_[istart]), oversize2[0], MPI_DOUBLE, neighbor_[0][iNeighbor], 0, SMILEI_COMM_1D, &(srequest[iNeighbor]) );
			//cout << "SUM : " << smilei_rk << " send " << oversize2[0] << " data to " << neighbor_[0][iNeighbor] << " starting at " << istart << endl;
		} // END of Send

		if (neighbor_[0][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
			istart = ( (iNeighbor+1)%2 ) * ( n_elem[0]- oversize2[0] ) + (1-(iNeighbor+1)%2) * ( 0 );
			MPI_Irecv( &( (buf[(iNeighbor+1)%2]).data_[0] ), oversize2[0], MPI_DOUBLE, neighbor_[0][(iNeighbor+1)%2], 0, SMILEI_COMM_1D, &(rrequest[(iNeighbor+1)%2]) );
			//cout << "SUM : " << smilei_rk << " recv " << oversize2[0] << " data to " << neighbor_[0][(iNeighbor+1)%2] << " starting at " << istart << endl;
		} // END of Recv

	} // END for iNeighbor


	for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
		if (neighbor_[0][iNeighbor]!=MPI_PROC_NULL ) {
			MPI_Wait( &(srequest[iNeighbor]), &(sstat[iNeighbor]) );
		}
		if (neighbor_[0][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
			MPI_Wait( &(rrequest[(iNeighbor+1)%2]), &(rstat[(iNeighbor+1)%2]) );
		}
	}


	// Synchro before summing, to not sum with data ever sum
	barrier();
	/********************************************************************************/
	// Sum data on each process, same operation on both side
	/********************************************************************************/
	for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
		istart = ( (iNeighbor+1)%2 ) * ( n_elem[0]- oversize2[0] ) + (1-(iNeighbor+1)%2) * ( 0 );
		// Using Receiver point of vue
		if (neighbor_[0][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
		  //cout << "SUM : " << smilei_rk << " sum " << oversize2[0] << " data from " << istart << endl;
			for (unsigned int i=0 ; i<oversize2[0] ; i++)
				f1D->data_[istart+i] += (buf[(iNeighbor+1)%2])(i);
		}
	} // END for iNeighbor


} // END sumField


void SmileiMPI_Cart1D::exchangeField( Field* field )
{
	std::vector<unsigned int> n_elem   = field->dims_;
	std::vector<unsigned int> isPrimal = field->isPrimal_;
	Field1D* f1D =  static_cast<Field1D*>(field);

	// Loop over dimField
	// See sumField for details
	int istart;
	MPI_Status sstat[2];
	MPI_Status rstat[2];
	MPI_Request srequest[2];
	MPI_Request rrequest[2];
	// Loop over neighbors in a direction
	// Send to neighbor_[0][iNeighbor] / Recv from neighbor_[0][(iNeighbor+1)%2] :
	// See in exchangeParticles()
	for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {

		if (neighbor_[0][iNeighbor]!=MPI_PROC_NULL) {
			istart = iNeighbor * ( n_elem[0]- (2*oversize[0]+1+isPrimal[0]) ) + (1-iNeighbor) * ( 2*oversize[0]+1-(1-isPrimal[0]) );
			MPI_Isend( &(f1D->data_[istart]), 1, MPI_DOUBLE, neighbor_[0][iNeighbor], 0, SMILEI_COMM_1D, &(srequest[iNeighbor]) );
			//cout << "EXCH : " << smilei_rk << " send " << oversize[0] << " data to " << neighbor_[0][iNeighbor] << " starting at " << istart << endl;
		} // END of Send

		if (neighbor_[0][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
			istart = ( (iNeighbor+1)%2 ) * ( n_elem[0] - 1 ) + (1-(iNeighbor+1)%2) * ( 0 )  ;
			MPI_Irecv( &(f1D->data_[istart]), 1, MPI_DOUBLE, neighbor_[0][(iNeighbor+1)%2], 0, SMILEI_COMM_1D, &(rrequest[(iNeighbor+1)%2]) );
			//cout << "EXCH : " << smilei_rk << " recv " << oversize[0] << " data to " << neighbor_[0][(iNeighbor+1)%2] << " starting at " << istart << endl;

		} // END of Recv

	} // END for iNeighbor


	for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
		if (neighbor_[0][iNeighbor]!=MPI_PROC_NULL) {
			MPI_Wait( &(srequest[iNeighbor]), &(sstat[iNeighbor]) );
		}
		if (neighbor_[0][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
			MPI_Wait( &(rrequest[(iNeighbor+1)%2]), &(rstat[(iNeighbor+1)%2]) );
		}
	}



} // END exchangeField


void SmileiMPI_Cart1D::writeField( Field* field, string name )
{
	Field1D* f1D =  static_cast<Field1D*>(field);
	std::vector<unsigned int> n_elem = field->dims_;
	int istart = oversize[0];
	if (smilei_rk!=0) istart+=1;  // f1D_current[n_elem[0]-2*oversize[0]+1] = f1D_west[oversize[0]]
	int bufsize = n_elem[0]- 2*oversize[0] - f1D->isPrimal_[0];
    
	if (smilei_rk!=0) {
		if (f1D->isPrimal_[0] == 0) bufsize-=1;
		else if (smilei_rk!=smilei_sz-1) bufsize-=1;
	}


	std::ofstream ff;

	for ( int i_rk = 0 ; i_rk < smilei_sz ; i_rk++ ) {
		if (i_rk==smilei_rk) {
			if (smilei_rk==0) ff.open(name.c_str(), ios::out);
			else ff.open(name.c_str(), ios::app);
			//cout << i_rk << " write " << bufsize-1 << " elements from " << istart << " to " << istart+bufsize-1 <<  endl;
			for (int i=istart ; i<istart+bufsize ; i++)
				ff << f1D->data_[i] << endl;
			if (smilei_rk==smilei_sz-1)ff << endl;
			//if (smilei_rk==smilei_sz-1)ff << endl;
			ff.close();
		}
		barrier();
	}


} // END writeField

