#include "SmileiMPI_Cart1D.h"

#include <cmath>
#include <cstring>

#include <string>

#include <mpi.h>
#include "Species.h"

#include "ElectroMagn.h"
#include "Field1D.h"

#include "Tools.h"

#include "Field2D.h"

using namespace std;

// ---------------------------------------------------------------------------------------------------------------------
// SmileiMPI_Cart1D: creator for Smilei MPI environment in 1D cartesian
// ---------------------------------------------------------------------------------------------------------------------
SmileiMPI_Cart1D::SmileiMPI_Cart1D( int* argc, char*** argv )
: SmileiMPI( argc, argv )
{
}


// ---------------------------------------------------------------------------------------------------------------------
// SmileiMPI_Cart1D: creator for Smilei MPI environment in 1D cartesian
// ---------------------------------------------------------------------------------------------------------------------
SmileiMPI_Cart1D::SmileiMPI_Cart1D( SmileiMPI* smpi)
: SmileiMPI( smpi )
{
    ndims_ = 1;
    number_of_procs  = new int[ndims_];
    coords_  = new int[ndims_];
    periods_  = new int[ndims_];
    reorder_ = 0;
    
    nbNeighbors_ = 2;
    
    for (int i=0 ; i<ndims_ ; i++) periods_[i] = 0;
    for (int i=0 ; i<ndims_ ; i++) coords_[i] = 0;
    for (int i=0 ; i<ndims_ ; i++) number_of_procs[i] = 1;
    
    for (int iDim=0 ; iDim<ndims_ ; iDim++) {
        for (int i=0 ; i<nbNeighbors_ ; i++) {
            neighbor_[iDim][i] = MPI_PROC_NULL;
            buff_index_send[i].resize(0);
            buff_index_recv_sz[i] = 0;
        }
    }
    
}


// ---------------------------------------------------------------------------------------------------------------------
// SmileiMPI_Cart1D: creator for Smilei MPI environment in 1D cartesian
// ---------------------------------------------------------------------------------------------------------------------
SmileiMPI_Cart1D::~SmileiMPI_Cart1D()
{
    delete [] number_of_procs;
    delete [] periods_;
    delete [] coords_;
    
    if ( SMILEI_COMM_1D != MPI_COMM_NULL) MPI_Comm_free(&SMILEI_COMM_1D);
    
}
// ---------------------------------------------------------------------------------------------------------------------
// SmileiMPI_Cart1D: create the topology for Smilei MPI environment in 1D cartesian
// ---------------------------------------------------------------------------------------------------------------------
void SmileiMPI_Cart1D::createTopology(Params& params)
{

    for (unsigned int i=0 ; i<params.nDim_field ; i++) {
        params.n_space_global[i] = round(params.sim_length[i]/params.cell_length[i]);
        MESSAGE("Total number of cells in direction " << i << ": " << params.n_space_global[i]);
    }
    
    number_of_procs[0] = smilei_sz;

    // Geometry periodic in x
    if ( (params.bc_em_type_x[0]=="periodic") || (params.bc_em_type_x[1]=="periodic") ) {
        periods_[0] = 1;
        MESSAGE("Periodic geometry in x-direction");
    }  

    MPI_Cart_create( SMILEI_COMM_WORLD, ndims_, number_of_procs, periods_, reorder_, &SMILEI_COMM_1D );
    MPI_Cart_coords( SMILEI_COMM_1D, smilei_rk, ndims_, coords_ );
    
    
    for (int iDim=0 ; iDim<ndims_ ; iDim++) {
        MPI_Cart_shift( SMILEI_COMM_1D, iDim, 1, &(neighbor_[iDim][0]), &(neighbor_[iDim][1]) );
        //DEBUG(3, smilei_rk, "Neighbors of process in direction " << iDim << " : " << neighbor_[iDim][0] << " ; " << neighbor_[iDim][1] << " Null :" << MPI_PROC_NULL );
    }
    
    
    for (unsigned int i=0 ; i<params.nDim_field ; i++) {
        
        n_space_global[i] = params.n_space_global[i];
        if ( (!params.nspace_win_x)||(i!=0) ) {
            
            params.n_space[i] = params.n_space_global[i] / number_of_procs[i];
            cell_starting_global_index[i] = coords_[i]*(params.n_space_global[i] / number_of_procs[i]);
            
            if ( number_of_procs[i]*params.n_space[i] != params.n_space_global[i] ) {
                // Correction on the last MPI process of the direction to use the wished number of cells
                if (coords_[i]==number_of_procs[i]-1) {
                    params.n_space[i] = params.n_space_global[i] - params.n_space[i]*(number_of_procs[i]-1);
                }
            }
        }
        else { // if use_moving_window
            // Number of space in window (not split)
            params.n_space[i] = params.nspace_win_x / number_of_procs[i];
            cell_starting_global_index[i] = coords_[i]*(params.nspace_win_x / number_of_procs[i]);
            
            if ( number_of_procs[i]*(int)params.n_space[i] != params.nspace_win_x ) {
                // Correction on the last MPI process of the direction to use the wished number of cells
                if (coords_[i]==number_of_procs[i]-1) {
                    params.n_space[i] = params.nspace_win_x - params.n_space[i]*(number_of_procs[i]-1);
                }
            }
        }
        
        oversize[i] = params.oversize[i] = params.interpolation_order + (params.exchange_particles_each-1);
        if ( params.n_space[i] <= 2*oversize[i] ) {
            WARNING ( "Increase space resolution or reduce number of MPI process in direction " << i );
        }
        
        // min/max_local : describe local domain in which particles cat be moved
        //                 different from domain on which E, B, J are defined
        min_local[i] = (cell_starting_global_index[i]                  )*params.cell_length[i];
        max_local[i] = (cell_starting_global_index[i]+params.n_space[i])*params.cell_length[i];
        
        cell_starting_global_index[i] -= params.oversize[i];
        
    }
    
    //DEBUG(3, smilei_rk, "n_space = " << params.n_space[0] );
    
    
    // -------------------------------------------------------
    // Compute & store the ranks of processes dealing with the
    // corner of the simulation box
    // -------------------------------------------------------
    
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

    
}

void SmileiMPI_Cart1D::exchangeParticles(Species* species, int ispec, Params& params,int tnum, int iDim)
{
    
    Particles &cuParticles = species->particles;
    std::vector<int>* cubmin = &species->bmin;
    std::vector<int>* cubmax = &species->bmax;
    
    MPI_Status Stat;
    int n_particles;
    int tid;
    int tmp = 0;
    int k=0;
    int i,ii, iPart;
    int n_part_recv, n_part_send;
    
    // ------------------------------------------------------------------------------
    // Build lists of indexes of particle to exchange per neighbor
    // Computed from indexes_of_particles_to_exchange computed during particles' BC
    // ------------------------------------------------------------------------------
    
    std::vector< std::vector<int> >* indexes_of_particles_to_exchange_per_thd = &species->indexes_of_particles_to_exchange_per_thd;
    //std::vector<int>                 indexes_of_particles_to_exchange;
    
#pragma omp single
    {
        indexes_of_particles_to_exchange.clear();
    }
#pragma omp barrier 
    
    for (tid=0 ; tid < tnum ; tid++){
        tmp += ((*indexes_of_particles_to_exchange_per_thd)[tid]).size(); //Compute the position where to start copying
    }

    if (tnum == (int)indexes_of_particles_to_exchange_per_thd->size()-1){ //If last thread
        indexes_of_particles_to_exchange.resize( tmp + ((*indexes_of_particles_to_exchange_per_thd)[tnum]).size());
    }
#pragma omp barrier 
    //Copy the list per_thread to the global list
    //One thread at a time (works)
#pragma omp master
    {
        for (tid=0 ; tid < (int)indexes_of_particles_to_exchange_per_thd->size() ; tid++) {
            memcpy(&indexes_of_particles_to_exchange[k], &((*indexes_of_particles_to_exchange_per_thd)[tid])[0],((*indexes_of_particles_to_exchange_per_thd)[tid]).size()*sizeof(int));
            k += ((*indexes_of_particles_to_exchange_per_thd)[tid]).size();   
        }
        // All threads together (doesn't work)
        /*if (((*indexes_of_particles_to_exchange_per_thd)[tnum]).size() > 0){
         //cout << "tmp = "<<tmp << endl;
         //cout << "tnum = "<< tnum << endl;
         memcpy(&indexes_of_particles_to_exchange[tmp], &((*indexes_of_particles_to_exchange_per_thd)[tnum])[0],((*indexes_of_particles_to_exchange_per_thd)[tnum]).size()*sizeof(int));
         }*/
        //#pragma omp master
        //{
        sort( indexes_of_particles_to_exchange.begin(), indexes_of_particles_to_exchange.end() );
        
        n_part_send = indexes_of_particles_to_exchange.size();
        
        
        for (i=0 ; i<n_part_send ; i++) {
            iPart = indexes_of_particles_to_exchange[i];
            if      ( cuParticles.position(0,iPart) < min_local[0]) {
                buff_index_send[0].push_back( indexes_of_particles_to_exchange[i] );
            }
            else if ( cuParticles.position(0,iPart) >= max_local[0]) {
                buff_index_send[1].push_back( indexes_of_particles_to_exchange[i] );
            }
        } // END for iPart = f(i)
        
        Particles partVectorSend[1][2];
        partVectorSend[0][0].initialize(0, params, ispec);
        partVectorSend[0][1].initialize(0, params, ispec);
        Particles partVectorRecv[1][2];
        partVectorRecv[0][0].initialize(0, params, ispec);
        partVectorRecv[0][1].initialize(0, params, ispec);
        
        /********************************************************************************/
        // Exchange particles
        /********************************************************************************/
        // Loop over neighbors in a direction
        
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
            n_part_send = buff_index_send[iNeighbor].size();
            if ( (neighbor_[0][0]!=MPI_PROC_NULL) && (neighbor_[0][1]!=MPI_PROC_NULL) ) {
                //Send-receive
                MPI_Sendrecv( &n_part_send, 1, MPI_INT, neighbor_[0][iNeighbor], 0, &buff_index_recv_sz[(iNeighbor+1)%2], 1, MPI_INT, neighbor_[0][(iNeighbor+1)%2], 0, SMILEI_COMM_1D,&Stat);
            } else if (neighbor_[0][iNeighbor]!=MPI_PROC_NULL) {
                //Send
                MPI_Send( &n_part_send, 1, MPI_INT, neighbor_[0][iNeighbor], 0, SMILEI_COMM_1D);
            } else if (neighbor_[0][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
                //Receive
                MPI_Recv( &buff_index_recv_sz[(iNeighbor+1)%2], 1, MPI_INT, neighbor_[0][(iNeighbor+1)%2], 0, SMILEI_COMM_1D, &Stat);
            }
        }
        
        /********************************************************************************/
        // Define buffers to exchange buff_index_send[iNeighbor].size();
        /********************************************************************************/
        //! \todo Define this as a main parameter for the code so that it needs not be defined all the time
        
        /********************************************************************************/
        // Proceed to effective Particles' communications
        /********************************************************************************/
        
        //Number of properties per particles = nDim_Particles + 3 + 1 + 1
        int nbrOfProp( 6 );
        MPI_Datatype typePartSend, typePartRecv;
        
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
            n_part_send = buff_index_send[iNeighbor].size();
            n_part_recv = buff_index_recv_sz[(iNeighbor+1)%2];
            if ( (neighbor_[0][0]!=MPI_PROC_NULL) && (neighbor_[0][1]!=MPI_PROC_NULL) && (n_part_send!=0) && (n_part_recv!=0) ) {
                //Send-receive
                double x_max = params.cell_length[0]*( params.n_space_global[0] );
                for (int iPart=0 ; iPart<n_part_send ; iPart++) {
                    // Enabled periodicity in X
                    if ( ( iNeighbor==0 ) &&  (coords_[0] == 0 ) &&( cuParticles.position(0,buff_index_send[iNeighbor][iPart]) < 0. ) ) {
                        cuParticles.position(0,buff_index_send[iNeighbor][iPart])     += x_max;
                    }
                    else if ( ( iNeighbor==1 ) &&  (coords_[0] == number_of_procs[0]-1 ) && ( cuParticles.position(0,buff_index_send[iNeighbor][iPart]) >= x_max ) ) {
                        cuParticles.position(0,buff_index_send[iNeighbor][iPart])     -= x_max;
                    }
                    cuParticles.cp_particle(buff_index_send[iNeighbor][iPart], partVectorSend[0][iNeighbor]);
                }
                
                typePartSend = createMPIparticles( &(partVectorSend[0][iNeighbor]), nbrOfProp );
                
                partVectorRecv[0][(iNeighbor+1)%2].initialize( n_part_recv, params, ispec );
                typePartRecv = createMPIparticles( &(partVectorRecv[0][(iNeighbor+1)%2]), nbrOfProp );
                
                MPI_Sendrecv(&((partVectorSend[0][iNeighbor      ]).position(0,0)),        1, typePartSend, neighbor_[0][iNeighbor      ], 0,
                             &((partVectorRecv[0][(iNeighbor+1)%2]).position(0,0)),        1, typePartRecv, neighbor_[0][(iNeighbor+1)%2], 0, SMILEI_COMM_1D, &Stat);
                MPI_Type_free( &typePartSend );
                MPI_Type_free( &typePartRecv );
                
            } else if ( (neighbor_[0][iNeighbor]!=MPI_PROC_NULL) && (n_part_send!=0) ) {
                //Send
                partVectorSend[0][iNeighbor].reserve(n_part_send, 1);
                double x_max = params.cell_length[0]*( params.n_space_global[0] );
                for (int iPart=0 ; iPart<n_part_send ; iPart++) {
                    // Enabled periodicity in X
                    if ( ( iNeighbor==0 ) &&  (coords_[0] == 0 ) &&( cuParticles.position(0,buff_index_send[iNeighbor][iPart]) < 0. ) ) {
                        cuParticles.position(0,buff_index_send[iNeighbor][iPart])     += x_max;
                    }
                    else if ( ( iNeighbor==1 ) &&  (coords_[0] == number_of_procs[0]-1 ) && ( cuParticles.position(0,buff_index_send[iNeighbor][iPart]) >= x_max ) ) {
                        cuParticles.position(0,buff_index_send[iNeighbor][iPart])     -= x_max;
                    }
                    cuParticles.cp_particle(buff_index_send[iNeighbor][iPart], partVectorSend[0][iNeighbor]);
                }
                typePartSend = createMPIparticles( &(partVectorSend[0][iNeighbor]), nbrOfProp );
                MPI_Send( &((partVectorSend[0][iNeighbor]).position(0,0)), 1, typePartSend, neighbor_[0][iNeighbor], 0, SMILEI_COMM_1D);
                MPI_Type_free( &typePartSend );
                
            } else if ( (neighbor_[0][(iNeighbor+1)%2]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
                //Receive
                partVectorRecv[0][(iNeighbor+1)%2].initialize( buff_index_recv_sz[(iNeighbor+1)%2], params, ispec );
                typePartRecv = createMPIparticles( &(partVectorRecv[0][(iNeighbor+1)%2]), nbrOfProp );
                MPI_Recv( &((partVectorRecv[0][(iNeighbor+1)%2]).position(0,0)), 1, typePartRecv,  neighbor_[0][(iNeighbor+1)%2], 0, SMILEI_COMM_1D, &Stat );
                MPI_Type_free( &typePartRecv );
                
            }
        }
        
        /********************************************************************************/
        // Delete Particles included in buff_send/buff_recv
        /********************************************************************************/
        // Push lost particles at the end of bins
        //! \todo For loop on bins, can use openMP here.
        for (unsigned int ibin = 0 ; ibin < (*cubmax).size() ; ibin++ ) {
            //        DEBUG(ibin << " bounds " << (*cubmin)[ibin] << " " << (*cubmax)[ibin]);
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
                    cuParticles.overwrite_part1D((*cubmax)[ibin]-1, iPart );
                    (*cubmax)[ibin]--;
                    ii--;
                    iPart = indexes_of_particles_to_exchange[ii];
                }
                if (iPart >= (*cubmin)[ibin] && iPart < (*cubmax)[ibin]) { //On traite la dernière particule (qui peut aussi etre la premiere)
                    cuParticles.overwrite_part1D((*cubmax)[ibin]-1, iPart );
                    (*cubmax)[ibin]--;
                }
            }
        }
        //Shift the bins in memory
        //Warning: this loop must be executed sequentially. Do not use openMP here.
        for (int unsigned ibin = 1 ; ibin < (*cubmax).size() ; ibin++ ) { //First bin don't need to be shifted
            ii = (*cubmin)[ibin]-(*cubmax)[ibin-1]; // Shift the bin in memory by ii slots.
            iPart = min(ii,(*cubmax)[ibin]-(*cubmin)[ibin]); // Number of particles we have to shift = min (Nshift, Nparticle in the bin)
            if(iPart > 0) cuParticles.overwrite_part1D((*cubmax)[ibin]-iPart,(*cubmax)[ibin-1],iPart);
            (*cubmax)[ibin] -= ii;
            (*cubmin)[ibin] = (*cubmax)[ibin-1];
        }
        
        
        // Delete useless Particles
        //Theoretically, not even necessary to do anything as long you use bmax as the end of your iterator on particles.
        //Nevertheless, you might want to free memory and have the actual number of particles
        //really equal to the size of the vector. So we do:
        cuParticles.erase_particle_trail((*cubmax).back());
        //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        
        
        /********************************************************************************/
        // Clean lists of indexes of particle to exchange per neighbor
        /********************************************************************************/
        for (int i=0 ; i<nbNeighbors_ ; i++)
            buff_index_send[i].clear();
        /********************************************************************************/
        // Copy newly arrived particles back to the vector
        /********************************************************************************/
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
            
            n_part_recv = buff_index_recv_sz[(iNeighbor+1)%2];
            if ( (neighbor_[0][(iNeighbor+1)%2]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
                if (iNeighbor == 0) { // Copy particles coming from the right at the end of Particles Array
                    n_particles = species->getNbrOfParticles();
                    partVectorRecv[0][(iNeighbor+1)%2].cp_particles(n_part_recv, cuParticles,n_particles);
                    (*cubmax)[(*cubmax).size()-1] += n_part_recv ;
                } else {// Copy particles coming from the left at the beginning of Particles Array
                    //New particles are inserted at the end of bin 0 instead of begining to minimize data movement.
                    partVectorRecv[0][(iNeighbor+1)%2].cp_particles(n_part_recv, cuParticles,(*cubmax)[0]);
                    (*cubmax)[0] += n_part_recv ;
                    for (unsigned int ibin=1 ; ibin < (*cubmax).size() ; ibin++ ) {
                        (*cubmax)[ibin] += n_part_recv ;
                        (*cubmin)[ibin] = (*cubmax)[ibin-1] ;
                    }
                }
            }
        }
    } // END omp master 
    //DEBUG( 2, "\tProcess " << smilei_rk << " : " << species->getNbrOfParticles() << " Particles of species " << ispec );
} // END exchangeParticles


MPI_Datatype SmileiMPI_Cart1D::createMPIparticles( Particles* particles, int nbrOfProp )
{
    MPI_Datatype typeParticlesMPI;


    int nbrOfProp2(nbrOfProp);
    if (particles->isTestParticles) nbrOfProp2++;

    MPI_Aint address[nbrOfProp2];
    MPI_Get_address( &(particles->position(0,0)), &(address[0]) );
    MPI_Get_address( &(particles->momentum(0,0)), &(address[1]) );
    MPI_Get_address( &(particles->momentum(1,0)), &(address[2]) );
    MPI_Get_address( &(particles->momentum(2,0)), &(address[3]) );
    MPI_Get_address( &(particles->weight(0)),     &(address[4]) );
    MPI_Get_address( &(particles->charge(0)),     &(address[5]) );
    //MPI_Get_address( &(particles.position_old(0,0)), &address[6] )
    if (particles->isTestParticles)
        MPI_Get_address( &(particles->id(0)),     &(address[nbrOfProp2-1]) );

    int nbr_parts[nbrOfProp2];
    MPI_Aint disp[nbrOfProp2];
    MPI_Datatype partDataType[nbrOfProp2];

    for (int i=0 ; i<nbrOfProp2 ; i++)
        nbr_parts[i] = particles->size();
    disp[0] = 0;
    for (int i=1 ; i<nbrOfProp2 ; i++)
        disp[i] = address[i] - address[0];
    for (int i=0 ; i<nbrOfProp2 ; i++)
        partDataType[i] = MPI_DOUBLE;
    partDataType[nbrOfProp-1] = MPI_SHORT;
    if (particles->isTestParticles)
        partDataType[nbrOfProp2-1] = MPI_UNSIGNED;

    MPI_Type_struct( nbrOfProp2, &(nbr_parts[0]), &(disp[0]), &(partDataType[0]), &typeParticlesMPI);
    MPI_Type_commit( &typeParticlesMPI );

    return typeParticlesMPI;
} // END createMPIparticles


void SmileiMPI_Cart1D::sumField( Field* field )
{
    std::vector<unsigned int> n_elem = field->dims_;
    Field1D* f1D =  static_cast<Field1D*>(field);
    
    // Use a buffer per direction to exchange data before summing
    Field1D buf[ nbNeighbors_ ];
    // Size buffer is 2 oversize (1 inside & 1 outside of the current subdomain)
    std::vector<unsigned int> oversize2 = oversize;
    oversize2[0] *= 2;
                
    oversize2[0] += 1 + f1D->isDual_[0];
    for (int i=0; i<nbNeighbors_ ; i++)  buf[i].allocateDims( oversize2 );
    
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
    std::vector<unsigned int> isDual = field->isDual_;
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
            istart = iNeighbor * ( n_elem[0]- (2*oversize[0]+1+isDual[0]) ) + (1-iNeighbor) * ( 2*oversize[0]+1-(1-isDual[0]) );
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
void SmileiMPI_Cart1D::exchangeField_movewin( Field* field, int clrw )
{
    std::vector<unsigned int> n_elem   = field->dims_;
    std::vector<unsigned int> isDual = field->isDual_;
    Field1D* f1D =  static_cast<Field1D*>(field);
    int istart, /*iDim,*/ iNeighbor, bufsize;
    void* b;
    
//    iDim = 0; // We exchange only in the X direction for movewin.
    iNeighbor = 0; // We send only towards the West and receive from the East.

    bufsize = clrw * sizeof(double) + 2 * MPI_BSEND_OVERHEAD; //Max number of doubles in the buffer. Careful, there might be MPI overhead to take into account.
    b=(void *)malloc(bufsize);
    MPI_Buffer_attach( b, bufsize);        
    
    // Loop over dimField
    // See sumField for details
    MPI_Status rstat;
    MPI_Request rrequest;
        

    if (neighbor_[0][iNeighbor]!=MPI_PROC_NULL) {
        istart = 2*oversize[0] + 1 + isDual[0]  ;
        MPI_Bsend( &(f1D->data_[istart]), clrw, MPI_DOUBLE, neighbor_[0][iNeighbor], 0, SMILEI_COMM_1D);
    } // END of Send

    field->shift_x(clrw);
        
    if (neighbor_[0][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
        istart = ( (iNeighbor+1)%2 ) * ( n_elem[0] - clrw ) + (1-(iNeighbor+1)%2) * ( 0 )  ;
        MPI_Irecv( &(f1D->data_[istart]), clrw, MPI_DOUBLE, neighbor_[0][(iNeighbor+1)%2], 0, SMILEI_COMM_1D, &rrequest );
    } // END of Recv
    
    
    if (neighbor_[0][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
        MPI_Wait( &rrequest, &rstat );
    }
    
    MPI_Buffer_detach( &b, &bufsize);        
    free(b);
    
    
} // END exchangeField
