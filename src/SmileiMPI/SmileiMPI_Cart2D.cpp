#include "SmileiMPI_Cart2D.h"

#include <cmath>
#include <string>

#include <mpi.h>

#include "Species.h"

#include "ElectroMagn.h"
#include "Field2D.h"

#include "Tools.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// SmileiMPI_Cart2D: creator for Smilei MPI environment in 2D cartesian geometry
// ---------------------------------------------------------------------------------------------------------------------
SmileiMPI_Cart2D::SmileiMPI_Cart2D( int* argc, char*** argv )
: SmileiMPI( argc, argv )
{
}


// ---------------------------------------------------------------------------------------------------------------------
// SmileiMPI_Cart2D: creator for Smilei MPI environment in 2D cartesian geometry
// ---------------------------------------------------------------------------------------------------------------------
SmileiMPI_Cart2D::SmileiMPI_Cart2D( SmileiMPI* smpi)
: SmileiMPI( smpi )
{
    ndims_ = 2;
    number_of_procs = new int(ndims_);
    coords_  = new int(ndims_);
    periods_  = new int(ndims_);
    reorder_ = 0;
    
    nbNeighbors_ = 2; // number of neighbor processes per direction
    
    for (int i=0 ; i<ndims_ ; i++) periods_[i] = 0;
    for (int i=0 ; i<ndims_ ; i++) coords_[i] = 0;
    for (int i=0 ; i<ndims_ ; i++) number_of_procs[i] = 1;
    
    for (int iDim=0 ; iDim<ndims_ ; iDim++)
        for (int iNeighbors=0 ; iNeighbors<nbNeighbors_ ; iNeighbors++)
            neighbor_[iDim][iNeighbors] = MPI_PROC_NULL;
    
    for (int iDim=0 ; iDim<ndims_ ; iDim++) {
        for (int i=0 ; i<nbNeighbors_ ; i++) {
            buff_index_send[iDim][i].resize(0);
            buff_index_recv_sz[iDim][i] = 0;
        }
    }
    
}


// ---------------------------------------------------------------------------------------------------------------------
// SmileiMPI_Cart2D: creator for Smilei MPI environment in 2D cartesian geometry
// ---------------------------------------------------------------------------------------------------------------------
SmileiMPI_Cart2D::~SmileiMPI_Cart2D()
{
    for (int ix_isPrim=0 ; ix_isPrim<1 ; ix_isPrim++) {
        for (int iy_isPrim=0 ; iy_isPrim<1 ; iy_isPrim++) {
            MPI_Type_free( &ntype_[0][ix_isPrim][iy_isPrim]); //line
            MPI_Type_free( &ntype_[1][ix_isPrim][iy_isPrim]); // column
            MPI_Type_free( &ntype_[2][ix_isPrim][iy_isPrim]); // several lines
            
            MPI_Type_free( &ntypeSum_[0][ix_isPrim][iy_isPrim]); //line
            MPI_Type_free( &ntypeSum_[1][ix_isPrim][iy_isPrim]); // column
        }
    }
    
    delete number_of_procs;
    delete periods_;
    delete coords_;
    

    if ( SMILEI_COMM_2D != MPI_COMM_NULL) MPI_Comm_free(&SMILEI_COMM_2D);
}
// ---------------------------------------------------------------------------------------------------------------------
// SmileiMPI_Cart2D: create the topology for Smilei MPI environment in 2D cartesian geometry
// ---------------------------------------------------------------------------------------------------------------------
void SmileiMPI_Cart2D::createTopology(PicParams& params)
{
    for (unsigned int i=0 ; i<params.nDim_field ; i++) {
        params.n_space_global[i] = round(params.sim_length[i]/params.cell_length[i]);
        MESSAGE("Total number of cells in direction " << i << ": " << params.n_space_global[i]);
    }
    
    if (params.number_of_procs[0]!=0) {
        for (unsigned int i=0 ; i<params.nDim_field ; i++)
            number_of_procs[i] = params.number_of_procs[i];
        if (number_of_procs[0]*number_of_procs[1]!=smilei_sz) {
            DEBUG(3,"Domain decomposition specified in the namelist don't match with the number of MPI process");
            DEBUG(3,"\tit will be computed to be as square as possible");
            for (unsigned int i=0 ; i<params.nDim_field ; i++)
                params.number_of_procs[i] = 0;
        }
    }
    if (params.number_of_procs[0]==0) {
	double tmp(0.);
        //      without moving-window
        if (!params.nspace_win_x)
            tmp  = params.res_space[0]*params.sim_length[0] / ( params.res_space[1]*params.sim_length[1] );
        else // with moving-window
            tmp = params.nspace_win_x * 2.0 * M_PI / ( params.res_space[1]*params.sim_length[1] );

        number_of_procs[0] = min( smilei_sz, max(1, (int)sqrt ( (double)smilei_sz*tmp*tmp) ) );
        number_of_procs[1] = (int)(smilei_sz / number_of_procs[0]);
        
        while ( number_of_procs[0]*number_of_procs[1] != smilei_sz ) {
            if (number_of_procs[0]>=number_of_procs[1] ) {
                number_of_procs[0]++;
                number_of_procs[1] = (int)(smilei_sz / number_of_procs[0]);
            }
            else {
                number_of_procs[1]++;
                number_of_procs[0] = (int)(smilei_sz / number_of_procs[1]);
            }
        }
        
    }
    // Force configuration of MPI domain decomposition
    //number_of_procs[0] = 1;
    //number_of_procs[1] = 16;
    MESSAGE("MPI Domain decomposition : " << smilei_sz << " = " << number_of_procs[0] << " x " << number_of_procs[1]);
    MESSAGE("Boundary conditions in x- & y-directions: "<< params.bc_em_type_long << ", " << params.bc_em_type_trans);
            
    // Geometry periodic in x
    if (params.bc_em_type_long=="periodic") {
        periods_[0] = 1;
        MESSAGE(1,"applied topology for periodic BCs in x-direction");
    }
    // Geometry periodic in y
    if (params.bc_em_type_trans=="periodic") {
        periods_[1] = 1;
        MESSAGE(2,"applied topology for periodic BCs in y-direction");
    }
    MPI_Cart_create( SMILEI_COMM_WORLD, ndims_, number_of_procs, periods_, reorder_, &SMILEI_COMM_2D );
    MPI_Cart_coords( SMILEI_COMM_2D, smilei_rk, ndims_, coords_ );
    
    
    
    
    for (int iDim=0 ; iDim<ndims_ ; iDim++) {
        MPI_Cart_shift( SMILEI_COMM_2D, iDim, 1, &(neighbor_[iDim][0]), &(neighbor_[iDim][1]) );
        //DEBUG(3,smilei_rk,"Neighbors of process in direction " << iDim << " : " << neighbor_[iDim][0] << " - " << neighbor_[iDim][1]  );
    }
    
    
    for (unsigned int i=0 ; i<params.nDim_field ; i++) {
        
        n_space_global[i] = params.n_space_global[i];
        
        if ((!params.nspace_win_x)||(i!=0)) {
            
            params.n_space[i] = params.n_space_global[i] / number_of_procs[i];
            cell_starting_global_index[i] = coords_[i]*(params.n_space_global[i] / number_of_procs[i]);
            
            if ( number_of_procs[i]*params.n_space[i] != params.n_space_global[i] ) {
                // Correction on the last MPI process of the direction to use the wished number of cells
                if (coords_[i]==number_of_procs[i]-1) {
                    params.n_space[i] = params.n_space_global[i] - params.n_space[i]*(number_of_procs[i]-1);
                }
            }
        }
        else {  // if use_moving_window
            // Number of space in window (not split)
            params.n_space[i] =  params.nspace_win_x / number_of_procs[i];
            cell_starting_global_index[i] = coords_[i]*(params.nspace_win_x / number_of_procs[i]);
            
            if ( number_of_procs[i]*params.n_space[i] != params.nspace_win_x ) {
                // Correction on the last MPI process of the direction to use the wished number of cells
                if (coords_[i]==number_of_procs[i]-1) {
                    params.n_space[i] = params.nspace_win_x - params.n_space[i]*(number_of_procs[i]-1);
                }
            }
        }
        
        oversize[i] = params.oversize[i] = params.interpolation_order + (params.exchange_particles_each-1);
        if ( params.n_space[i] <= 2*oversize[i] ) {
            WARNING ( "Increase space resolution or reduce number of MPI process in direction " << i << " "<< params.n_space[i]); 
        }
        
        // min/max_local : describe local domain in which particles cat be moved
        //                 different from domain on which E, B, J are defined
        min_local[i] = (cell_starting_global_index[i]                  )*params.cell_length[i];
        max_local[i] = (cell_starting_global_index[i]+params.n_space[i])*params.cell_length[i];
        //PMESSAGE( 0, smilei_rk, "min_local / mac_local on " << smilei_rk << " = " << min_local[i] << " / " << max_local[i] << " selon la direction " << i );
        
        cell_starting_global_index[i] -= params.oversize[i];
        
    }
    
    //DEBUG(10,"n_space / rank " << smilei_rk << " = " << params.n_space[0] << " " << params.n_space[1] );
    
    extrem_ranks[0][0] = 0;
    int rank_min =  0;
    if ( (coords_[0] == 0) && (coords_[1] == 0) )
        rank_min = smilei_rk;
    MPI_Allreduce(&rank_min, &extrem_ranks[0][0], 1, MPI_INT, MPI_SUM, SMILEI_COMM_2D);
    
    extrem_ranks[0][1] = 0;
    int rank_max = 0;
    if ( (coords_[0]==0) && (coords_[1]==number_of_procs[1]-1) )
        rank_max = smilei_rk;
    MPI_Allreduce(&rank_max, &extrem_ranks[0][1], 1, MPI_INT, MPI_SUM, SMILEI_COMM_2D);
    
    extrem_ranks[1][0] = 0;
    rank_max = 0;
    if ( (coords_[1]==0) && (coords_[0]==number_of_procs[0]-1) )
        rank_max = smilei_rk;
    MPI_Allreduce(&rank_max, &extrem_ranks[1][0], 1, MPI_INT, MPI_SUM, SMILEI_COMM_2D);
    
    extrem_ranks[1][1] = 0;
    rank_max = 0;
    if ( (coords_[0]==number_of_procs[0]-1) && (coords_[1]==number_of_procs[1]-1) )
        rank_max = smilei_rk;
    MPI_Allreduce(&rank_max, &extrem_ranks[1][1], 1, MPI_INT, MPI_SUM, SMILEI_COMM_2D);
    
    
}

void SmileiMPI_Cart2D::exchangeParticles(Species* species, int ispec, PicParams& params, int tnum)
{
    Particles &cuParticles = species->particles;
    std::vector<int>* cubmin = &species->bmin;
    std::vector<int>* cubmax = &species->bmax;
    
    std::vector< std::vector<int> >* indexes_of_particles_to_exchange_per_thd = &species->indexes_of_particles_to_exchange_per_thd;
    std::vector<int>                 indexes_of_particles_to_exchange;
    
#pragma omp master
    {
        /********************************************************************************/
        // Build lists of indexes of particle to exchange per neighbor
        // Computed from indexes_of_particles_to_exchange computed during particles' BC
        /********************************************************************************/
        indexes_of_particles_to_exchange.clear();
        
        int tmp = 0;
        for (int tid=0 ; tid < indexes_of_particles_to_exchange_per_thd->size() ; tid++)
            tmp += ((*indexes_of_particles_to_exchange_per_thd)[tid]).size();
        indexes_of_particles_to_exchange.resize( tmp );
        
        int k=0;
        for (int tid=0 ; tid < indexes_of_particles_to_exchange_per_thd->size() ; tid++) {
            for (int ipart = 0 ; ipart < ((*indexes_of_particles_to_exchange_per_thd)[tid]).size() ; ipart++ ) {
                indexes_of_particles_to_exchange[k] =  (*indexes_of_particles_to_exchange_per_thd)[tid][ipart] ;
                k++;
            }
            ((*indexes_of_particles_to_exchange_per_thd))[tid].clear();
        }
        sort( indexes_of_particles_to_exchange.begin(), indexes_of_particles_to_exchange.end() );
        
        int n_part_send = indexes_of_particles_to_exchange.size();
        int n_part_recv;
        
        int ii,iPart;
        int n_particles,nmove,lmove;
        int shift[(*cubmax).size()+1];//how much we need to shift each bin in order to leave room for the new particles
        double dbin;
        
        dbin = params.cell_length[0]*params.clrw; //width of a bin.
        for (unsigned int j=0; j<(*cubmax).size()+1 ;j++){
            shift[j]=0;
        }
        
        for (int i=0 ; i<n_part_send ; i++) {
            iPart = indexes_of_particles_to_exchange[i];
            // Y direction managed firstly to force considering periodic modification
            for (int iDim=1 ; iDim>=0 ; iDim--) {
                if ( cuParticles.position(iDim,iPart) < min_local[iDim]) {
                    buff_index_send[iDim][0].push_back( indexes_of_particles_to_exchange[i] );
                    break;
                }
                if ( cuParticles.position(iDim,iPart) >= max_local[iDim]) {
                    buff_index_send[iDim][1].push_back( indexes_of_particles_to_exchange[i] );
                    break;
                }
            }
        } // END for iPart = f(i)
        
        Particles partVectorSend[2][2];
        partVectorSend[0][0].initialize(0,cuParticles.dimension());
        partVectorSend[0][1].initialize(0,cuParticles.dimension());
        partVectorSend[1][0].initialize(0,cuParticles.dimension());
        partVectorSend[1][1].initialize(0,cuParticles.dimension());
        Particles partVectorRecv[2][2];
        partVectorRecv[0][0].initialize(0,cuParticles.dimension());
        partVectorRecv[0][1].initialize(0,cuParticles.dimension());
        partVectorRecv[1][0].initialize(0,cuParticles.dimension());
        partVectorRecv[1][1].initialize(0,cuParticles.dimension());
        
        Particles diagonalParticles;
        diagonalParticles.initialize(0,cuParticles.dimension());
        
        /********************************************************************************/
        // Exchange particles
        /********************************************************************************/
        for (int iDim=0 ; iDim<ndims_ ; iDim++) {
            
            MPI_Status sstat    [ndims_][2];
            MPI_Status rstat    [ndims_][2];
            MPI_Request srequest[ndims_][2];
            MPI_Request rrequest[ndims_][2];
            
            /********************************************************************************/
            // Exchange number of particles to exchange to establish or not a communication
            /********************************************************************************/
            for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
                if (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {
                    n_part_send = (buff_index_send[iDim][iNeighbor]).size();
                    MPI_Isend( &n_part_send, 1, MPI_INT, neighbor_[iDim][iNeighbor], 0, SMILEI_COMM_2D, &(srequest[iDim][iNeighbor]) );
                } // END of Send
                if (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
                    MPI_Irecv( &(buff_index_recv_sz[iDim][(iNeighbor+1)%2]), 1, MPI_INT, neighbor_[iDim][(iNeighbor+1)%2], 0, SMILEI_COMM_2D, &(rrequest[iDim][(iNeighbor+1)%2]) );
                }
            }
            barrier();
            
            /********************************************************************************/
            // Wait for end of communications over number of particles
            /********************************************************************************/
            for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
                if (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {
                    MPI_Wait( &(srequest[iDim][iNeighbor]), &(sstat[iDim][iNeighbor]) );
                }
                if (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
                    MPI_Wait( &(rrequest[iDim][(iNeighbor+1)%2]), &(rstat[iDim][(iNeighbor+1)%2]) );
                    if (buff_index_recv_sz[iDim][(iNeighbor+1)%2]!=0) {
                        partVectorRecv[iDim][(iNeighbor+1)%2].initialize( buff_index_recv_sz[iDim][(iNeighbor+1)%2], cuParticles.dimension());
                    }
                }
            }
            barrier();
            
            /********************************************************************************/
            // Define buffers to exchange buff_index_send[iDim][iNeighbor].size();
            /********************************************************************************/
            
            
            /********************************************************************************/
            // Proceed to effective Particles' communications
            /********************************************************************************/

            // Number of properties per particles = nDim_Particles + 3 + 1 + 1
	    int nbrOfProp( 7 );
	    MPI_Datatype typePartSend, typePartRecv;

            for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
                
                // n_part_send : number of particles to send to current neighbor
                n_part_send = (buff_index_send[iDim][iNeighbor]).size();
                if ( (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) && (n_part_send!=0) ) {
                    double x_max = params.cell_length[0]*( params.n_space_global[0] );
                    double y_max = params.cell_length[1]*( params.n_space_global[1] );
                    for (int iPart=0 ; iPart<n_part_send ; iPart++) {

			if (iDim==0) {
			    // Enabled periodicity in X
			    if ( ( iNeighbor==0 ) &&  (coords_[0] == 0 ) &&( cuParticles.position(0,buff_index_send[iDim][iNeighbor][iPart]) < 0. ) ) {
                                cuParticles.position(0,buff_index_send[iDim][iNeighbor][iPart])     += x_max;
                            }
                            else if ( ( iNeighbor==1 ) &&  (coords_[0] == number_of_procs[0]-1 ) && ( cuParticles.position(0,buff_index_send[iDim][iNeighbor][iPart]) >= x_max ) ) {
                                cuParticles.position(0,buff_index_send[iDim][iNeighbor][iPart])     -= x_max;
                            }
			}
			else if (iDim==1) {
			    // Enabled periodicity in Y
                            if ( ( iNeighbor==0 ) &&  (coords_[1] == 0 ) &&( cuParticles.position(1,buff_index_send[iDim][iNeighbor][iPart]) < 0. ) ) {
                                cuParticles.position(1,buff_index_send[iDim][iNeighbor][iPart])     += y_max;
                            }
                            else if ( ( iNeighbor==1 ) &&  (coords_[1] == number_of_procs[1]-1 ) && ( cuParticles.position(1,buff_index_send[iDim][iNeighbor][iPart]) >= y_max ) ) {
                                cuParticles.position(1,buff_index_send[iDim][iNeighbor][iPart])     -= y_max;
                            }
			}
                        cuParticles.cp_particle(buff_index_send[iDim][iNeighbor][iPart], partVectorSend[iDim][iNeighbor]);
                    }

		    typePartSend = createMPIparticles( &(partVectorSend[iDim][iNeighbor]), nbrOfProp );
		    MPI_Isend( &((partVectorSend[iDim][iNeighbor]).position(0,0)), 1, typePartSend, neighbor_[iDim][iNeighbor], 0, SMILEI_COMM_2D, &(srequest[iDim][iNeighbor]) );
		    MPI_Type_free( &typePartSend );

                } // END of Send
                
                n_part_recv = buff_index_recv_sz[iDim][(iNeighbor+1)%2];
                if ( (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
		    typePartRecv = createMPIparticles( &(partVectorRecv[iDim][(iNeighbor+1)%2]), nbrOfProp );
		    MPI_Irecv( &((partVectorRecv[iDim][(iNeighbor+1)%2]).position(0,0)), 1, typePartRecv,  neighbor_[iDim][(iNeighbor+1)%2], 0, SMILEI_COMM_2D, &(rrequest[iDim][(iNeighbor+1)%2]) );
		    MPI_Type_free( &typePartRecv );

                } // END of Recv
                
            } // END for iNeighbor
            
            
            /********************************************************************************/
            // Wait for end of communications over Particles
            /********************************************************************************/
            for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
                
                n_part_send = buff_index_send[iDim][iNeighbor].size();
                n_part_recv = buff_index_recv_sz[iDim][(iNeighbor+1)%2];
                
                if ( (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) && (n_part_send!=0) ) {
                    MPI_Wait( &(srequest[iDim][iNeighbor]), &(sstat[iDim][iNeighbor]) );
                }
                
                if ( (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
                    MPI_Wait( &(rrequest[iDim][(iNeighbor+1)%2]), &(rstat[iDim][(iNeighbor+1)%2]) );
                    
                    // Extract corner particles, not managed in the following process 
                    // but reinjected at the end in the main list
                    for (int iPart=0 ; iPart<n_part_recv; iPart++ )
                        if ( !(partVectorRecv[iDim][(iNeighbor+1)%2]).is_part_in_domain(iPart, this) )
                            (partVectorRecv[iDim][(iNeighbor+1)%2]).cp_particle(iPart, diagonalParticles);
                    for (int iPart=n_part_recv-1 ; iPart>=0; iPart-- ) {
                        if ( !(partVectorRecv[iDim][(iNeighbor+1)%2]).is_part_in_domain(iPart, this) ) {
                            (partVectorRecv[iDim][(iNeighbor+1)%2]).erase_particle(iPart);
                            buff_index_recv_sz[iDim][(iNeighbor+1)%2]--;
                        }
                    }
                    
                }
                
            }
            barrier();
            /********************************************************************************/
            // Clean lists of indexes of particle to exchange per neighbor
            /********************************************************************************/
            for (int i=0 ; i<nbNeighbors_ ; i++)
                buff_index_send[iDim][i].clear();
            
        } //End of iDim loop
        
        /********************************************************************************/
        // Delete Particles included in buff_send/buff_recv
        /********************************************************************************/
        
        // Push lost particles at the end of bins
        //! \todo For loop on bins, can use openMP here.
        for (unsigned int ibin = 0 ; ibin < (*cubmax).size() ; ibin++ ) {
            //cout << "bounds " << (*cubmin)[ibin] << " " << (*cubmax)[ibin] << endl;
            ii = indexes_of_particles_to_exchange.size()-1;
            if (ii >= 0) { // Push lost particles to the end of the bin
                iPart = indexes_of_particles_to_exchange[ii];
                while (iPart >= (*cubmax)[ibin] && ii > 0) {
                    ii--;
                    iPart = indexes_of_particles_to_exchange[ii];
                }
                while (iPart == (*cubmax)[ibin]-1 && iPart >= (*cubmin)[ibin] && ii > 0) {
                    (*cubmax)[ibin]--;
                    ii--;
                    iPart = indexes_of_particles_to_exchange[ii];
                }
                while (iPart >= (*cubmin)[ibin] && ii > 0) {
                    cuParticles.overwrite_part2D((*cubmax)[ibin]-1, iPart );
                    (*cubmax)[ibin]--;
                    ii--;
                    iPart = indexes_of_particles_to_exchange[ii];
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
        // Delete useless Particles
        //Theoretically, not even necessary to do anything as long you use bmax as the end of your iterator on particles.
        //Nevertheless, you might want to free memory and have the actual number of particles
        //really equal to the size of the vector. So we do:
        cuParticles.erase_particle_trail((*cubmax).back());
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        //********************************************************************************/
        // Copy newly arrived particles back to the vector
        // WARNING: very different behaviour depending on which dimension particles are coming from.
        /********************************************************************************/
        //We first evaluate how many particles arrive in each bin. 
        //1) Count particles coming from south and north
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
            n_part_recv = buff_index_recv_sz[1][iNeighbor];
            for (unsigned int j=0; j<n_part_recv ;j++){
                ii = int((partVectorRecv[1][iNeighbor].position(0,j)-min_local[0])/dbin);//bin in which the particle goes.
                shift[ii+1]++; // It makes the next bins shift.
            }
        }
        //2) Add particles coming from west and east
        shift[1] += buff_index_recv_sz[0][0];//Particles coming from south all go to bin 0 and shift all the other bins.
        shift[(*cubmax).size()] += buff_index_recv_sz[0][1];//Used only to count the total number of particles arrived.
        
        //Evaluation of the necessary shift of all bins.
        //Must be done sequentially
        for (unsigned int j=1; j<(*cubmax).size()+1;j++){ //bin 0 is not shifted.Last element of shift stores total number of arriving particles.
            shift[j]+=shift[j-1];
        }
        //Make room for new particles
        //cuParticles.create_particles(shift[(*cubmax).size()]);
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
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
            n_part_recv = buff_index_recv_sz[0][iNeighbor];
            if ( (neighbor_[0][iNeighbor]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
                ii = iNeighbor*((*cubmax).size()-1);//0 if iNeighbor=0(particles coming from West) and (*cubmax).size()-1 otherwise.
                partVectorRecv[0][iNeighbor].overwrite_part2D(0, cuParticles,(*cubmax)[ii],n_part_recv);
                (*cubmax)[ii] += n_part_recv ;
            }
        }
        //iDim == 1) is the difficult case, when particles can arrive in any bin.
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
            n_part_recv = buff_index_recv_sz[1][iNeighbor];
            if ( (neighbor_[1][iNeighbor]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
                for(unsigned int j=0; j<n_part_recv; j++){
                    ii = int((partVectorRecv[1][iNeighbor].position(0,j)-min_local[0])/dbin);//bin in which the particle goes.
                    partVectorRecv[1][iNeighbor].overwrite_part2D(j, cuParticles,(*cubmax)[ii]);
                    (*cubmax)[ii] ++ ;
                }
            }
        }
        
        // Inject corner particles at the end of the list, update bmax
        for (int iPart = 0 ; iPart<diagonalParticles.size() ; iPart++) {
            diagonalParticles.cp_particle(iPart, cuParticles);
            (*indexes_of_particles_to_exchange_per_thd)[0].push_back(cuParticles.size()-1);
            (*cubmax)[(*cubmax).size()-1]++;
        }
        
    }//end of omp master
} // END exchangeParticles


MPI_Datatype SmileiMPI_Cart2D::createMPIparticles( Particles* particles, int nbrOfProp )
{
    MPI_Datatype typeParticlesMPI;

    MPI_Aint address[nbrOfProp];
    MPI_Get_address( &(particles->position(0,0)), &(address[0]) );
    MPI_Get_address( &(particles->position(1,0)), &(address[1]) );
    //MPI_Get_address( &(particles->position_old(0,0)), &(address[2]) );
    //MPI_Get_address( &(particles->position_old(1,0)), &(address[3]) );
    MPI_Get_address( &(particles->momentum(0,0)), &(address[2]) );
    MPI_Get_address( &(particles->momentum(1,0)), &(address[3]) );
    MPI_Get_address( &(particles->momentum(2,0)), &(address[4]) );
    MPI_Get_address( &(particles->weight(0)),     &(address[5]) );
    MPI_Get_address( &(particles->charge(0)),     &(address[6]) );

    int nbr_parts[nbrOfProp];
    MPI_Aint disp[nbrOfProp];
    MPI_Datatype partDataType[nbrOfProp];

    for (int i=0 ; i<nbrOfProp ; i++)
	nbr_parts[i] = particles->size();
    disp[0] = 0;
    for (int i=1 ; i<nbrOfProp ; i++)
	disp[i] = address[i] - address[0];
    for (int i=0 ; i<nbrOfProp ; i++)
	partDataType[i] = MPI_DOUBLE;
    partDataType[nbrOfProp-1] = MPI_SHORT;

    MPI_Type_struct( nbrOfProp, &(nbr_parts[0]), &(disp[0]), &(partDataType[0]), &typeParticlesMPI);
    MPI_Type_commit( &typeParticlesMPI );

    return typeParticlesMPI;
} // END createMPIparticles

void SmileiMPI_Cart2D::createType( PicParams& params )
{
    int nx0 = params.n_space[0] + 1 + 2*params.oversize[0];
    int ny0 = params.n_space[1] + 1 + 2*params.oversize[1];
    int clrw = params.clrw;
    
    // MPI_Datatype ntype_[nDim][primDual][primDual]
    int nx, ny;
    int nline, ncol;
    for (int ix_isPrim=0 ; ix_isPrim<2 ; ix_isPrim++) {
        nx = nx0 + ix_isPrim;
        for (int iy_isPrim=0 ; iy_isPrim<2 ; iy_isPrim++) {
            ny = ny0 + iy_isPrim;
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
            
        }
    }
    
} //END createType


void SmileiMPI_Cart2D::sumField( Field* field )
{
    std::vector<unsigned int> n_elem = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field2D* f2D =  static_cast<Field2D*>(field);
    
    
    // Use a buffer per direction to exchange data before summing
    Field2D buf[ndims_][ nbNeighbors_ ];
    // Size buffer is 2 oversize (1 inside & 1 outside of the current subdomain)
    std::vector<unsigned int> oversize2 = oversize;
    oversize2[0] *= 2;
    oversize2[0] += 1 + f2D->isDual_[0];
    oversize2[1] *= 2;
    oversize2[1] += 1 + f2D->isDual_[1];
    
    for (int iDim=0 ; iDim<ndims_ ; iDim++) {
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
            std::vector<unsigned int> tmp(ndims_,0);
            tmp[0] =    iDim  * n_elem[0] + (1-iDim) * oversize2[0];
            tmp[1] = (1-iDim) * n_elem[1] +    iDim  * oversize2[1];
            buf[iDim][iNeighbor].allocateDims( tmp );
        }
    }
    
    int istart, ix, iy;
    
    /********************************************************************************/
    // Send/Recv in a buffer data to sum
    /********************************************************************************/
    for (int iDim=0 ; iDim<ndims_ ; iDim++) {
        
        MPI_Datatype ntype = ntypeSum_[iDim][isDual[0]][isDual[1]];
        //      MPI_Status stat[2];
        //      MPI_Request request[2];
        MPI_Status sstat    [ndims_][2];
        MPI_Status rstat    [ndims_][2];
        MPI_Request srequest[ndims_][2];
        MPI_Request rrequest[ndims_][2];
        
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
            
            if (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {
                istart = iNeighbor * ( n_elem[iDim]- oversize2[iDim] ) + (1-iNeighbor) * ( 0 );
                ix = (1-iDim)*istart;
                iy =    iDim *istart;
                MPI_Isend( &(f2D->data_2D[ix][iy]), 1, ntype, neighbor_[iDim][iNeighbor], 0, SMILEI_COMM_2D, &(srequest[iDim][iNeighbor]) );
            } // END of Send
            
            if (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
                int tmp_elem = (buf[iDim][(iNeighbor+1)%2]).dims_[0]*(buf[iDim][(iNeighbor+1)%2]).dims_[1];
                MPI_Irecv( &( (buf[iDim][(iNeighbor+1)%2]).data_2D[0][0] ), tmp_elem, MPI_DOUBLE, neighbor_[iDim][(iNeighbor+1)%2], 0, SMILEI_COMM_2D, &(rrequest[iDim][(iNeighbor+1)%2]) );
            } // END of Recv
            
        } // END for iNeighbor
        
        
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
            if (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {
                MPI_Wait( &(srequest[iDim][iNeighbor]), &(sstat[iDim][iNeighbor]) );
            }
            if (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
                MPI_Wait( &(rrequest[iDim][(iNeighbor+1)%2]), &(rstat[iDim][(iNeighbor+1)%2]) );
            }
        }
        
        
        // Synchro before summing, to not sum with data ever sum
        // Merge loops, Sum direction by direction permits to not communicate with diagonal neighbors
        barrier();
        /********************************************************************************/
        // Sum data on each process, same operation on both side
        /********************************************************************************/
        
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
            istart = ( (iNeighbor+1)%2 ) * ( n_elem[iDim]- oversize2[iDim] ) + (1-(iNeighbor+1)%2) * ( 0 );
            int ix0 = (1-iDim)*istart;
            int iy0 =    iDim *istart;
            if (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
                for (unsigned int ix=0 ; ix< (buf[iDim][(iNeighbor+1)%2]).dims_[0] ; ix++) {
                    for (unsigned int iy=0 ; iy< (buf[iDim][(iNeighbor+1)%2]).dims_[1] ; iy++)
                        f2D->data_2D[ix0+ix][iy0+iy] += (buf[iDim][(iNeighbor+1)%2])(ix,iy);
                }
            } // END if
            
        } // END for iNeighbor
        
        barrier();
        
    } // END for iDim
    
} // END sumField


void SmileiMPI_Cart2D::exchangeField( Field* field )
{
    std::vector<unsigned int> n_elem   = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field2D* f2D =  static_cast<Field2D*>(field);
    
    int istart, ix, iy;
    
    // Loop over dimField
    for (int iDim=0 ; iDim<ndims_ ; iDim++) {
        
        MPI_Datatype ntype = ntype_[iDim][isDual[0]][isDual[1]];
        //      MPI_Status stat[2];
        //      MPI_Request request[2];
        MPI_Status sstat    [ndims_][2];
        MPI_Status rstat    [ndims_][2];
        MPI_Request srequest[ndims_][2];
        MPI_Request rrequest[ndims_][2];
        
        
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
            
            if (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {
                
                istart = iNeighbor * ( n_elem[iDim]- (2*oversize[iDim]+1+isDual[iDim]) ) + (1-iNeighbor) * ( 2*oversize[iDim] + isDual[iDim] );
                ix = (1-iDim)*istart;
                iy =    iDim *istart;
                MPI_Isend( &(f2D->data_2D[ix][iy]), 1, ntype, neighbor_[iDim][iNeighbor], 0, SMILEI_COMM_2D, &(srequest[iDim][iNeighbor]) );
                
            } // END of Send
            
            if (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
                
                istart = ( (iNeighbor+1)%2 ) * ( n_elem[iDim] - 1 ) + (1-(iNeighbor+1)%2) * ( 0 )  ;
                ix = (1-iDim)*istart;
                iy =    iDim *istart;
                MPI_Irecv( &(f2D->data_2D[ix][iy]), 1, ntype, neighbor_[iDim][(iNeighbor+1)%2], 0, SMILEI_COMM_2D, &(rrequest[iDim][(iNeighbor+1)%2]));
                
            } // END of Recv
            
        } // END for iNeighbor
        
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
            if (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {
                MPI_Wait( &(srequest[iDim][iNeighbor]), &(sstat[iDim][iNeighbor]) );
            }
            if (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
                MPI_Wait( &(rrequest[iDim][(iNeighbor+1)%2]), &(rstat[iDim][(iNeighbor+1)%2]) );
            }
        }
        
    } // END for iDim
    
    
} // END exchangeField

void SmileiMPI_Cart2D::exchangeField_movewin( Field* field, int clrw )
{
    std::vector<unsigned int> n_elem   = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field2D* f2D =  static_cast<Field2D*>(field);
    int istart, ix, iy, iDim, iNeighbor,bufsize;
    void* b;
    
    bufsize = clrw*n_elem[1]*sizeof(double)+ 2 * MPI_BSEND_OVERHEAD; //Max number of doubles in the buffer. Careful, there might be MPI overhead to take into account.
    b=(void *)malloc(bufsize);
    MPI_Buffer_attach( b, bufsize);        
    iDim = 0; // We exchange only in the X direction for movewin.
    iNeighbor = 0; // We send only towards the West and receive from the East.

    MPI_Datatype ntype = ntype_[2][isDual[0]][isDual[1]]; //ntype_[2] is clrw columns.
    MPI_Status rstat    ;
    MPI_Request rrequest;
    
    
    if (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {
        
        istart =  2*oversize[iDim] + 1 + isDual[iDim] ;
        ix = (1-iDim)*istart;
        iy =    iDim *istart;
        MPI_Bsend( &(f2D->data_2D[ix][iy]), 1, ntype, neighbor_[iDim][iNeighbor], 0, SMILEI_COMM_2D);
    } // END of Send
  
    //Once the message is in the buffer we can safely shift the field in memory. 
    field->shift_x(clrw);
    // and then receive the complementary field from the East.

    if (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
        
        istart = ( (iNeighbor+1)%2 ) * ( n_elem[iDim] - clrw ) + (1-(iNeighbor+1)%2) * ( 0 )  ;
        ix = (1-iDim)*istart;
        iy =    iDim *istart;
        MPI_Irecv( &(f2D->data_2D[ix][iy]), 1, ntype, neighbor_[iDim][(iNeighbor+1)%2], 0, SMILEI_COMM_2D, &rrequest);
    } // END of Recv
        
    
    if (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
        MPI_Wait( &rrequest, &rstat);
    }
    MPI_Buffer_detach( &b, &bufsize);        
    free(b);
    
    
} // END exchangeField_movewin
