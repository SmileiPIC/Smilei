#include "Patch.h"
#include <cstdlib>
#include <iostream>
#include <iomanip>

using namespace std;


Patch::Patch(PicParams& params, LaserParams& laser_params, SmileiMPI* smpi, unsigned int m0, unsigned int m1, unsigned int m2, unsigned int ipatch) {

        hindex = ipatch;
        if ( params.geometry == "1d3v" ) {
            mi.resize(1);
            Pcoordinates.resize(1);
            mi[0] = m0;
            Pcoordinates[0] = hindex;
        }
        else if ( params.geometry == "2d3v" ) {
            mi.resize(2);
            Pcoordinates.resize(2);
            mi[0] = m0;
            mi[1] = m1;
            compacthilbertindexinv(m0, m1, &Pcoordinates[0], &Pcoordinates[1], hindex);
        }
        else {
            mi.resize(3);
            Pcoordinates.resize(3);
            mi[0] = m0;
            mi[1] = m1;
            mi[2] = m2;
            compacthilbertindexinv(m0, m1, m2, &Pcoordinates[0], &Pcoordinates[1], &Pcoordinates[2], hindex);
        }

	Pcoordinates[1] = ipatch%m1;
	Pcoordinates[0] = ipatch/m1 ;
	//std::cout << "Coordonnées de " << ipatch << " : " << Pcoordinates[0] << " " << Pcoordinates[1] << std::endl;
	neighbor_.resize(params.nDim_field);
	for ( int iDim = 0 ; iDim < params.nDim_field ; iDim++ ) {
	    neighbor_[iDim].resize(2,-2);
	}
	if (Pcoordinates[0]>0)
	    neighbor_[0][0] = ipatch-m1;
	if (Pcoordinates[0]<m0-1)
	    neighbor_[0][1] = ipatch+m1;
	if (Pcoordinates[1]>0)
	    neighbor_[1][0] = ipatch-1;
	if (Pcoordinates[1]<m1-1)
	    neighbor_[1][1] = ipatch+1;

	// Manage y-periodicity only !
	SmileiMPI_Cart2D* smpi2D = static_cast<SmileiMPI_Cart2D*>(smpi);
	if ( (params.bc_em_type_trans=="periodic") ) {
	  if ( (smpi2D->getNbrOfProcs(1)==1) ) {
	    if ( (smpi2D->getProcCoord(1)==0) || (smpi2D->getProcCoord(1)==smpi2D->getNbrOfProcs(1)-1) ) {
	      if ( (Pcoordinates[1]==0) )
		neighbor_[1][0] = Pcoordinates[0]*m1+m1-1;
	      if ( (Pcoordinates[1]==m1-1) )
		neighbor_[1][1] = Pcoordinates[0]*m1+0;
	    }
	  }
	}
	
	//std::cout << "Voisin dir 0 : " << ipatch << " : " <<  neighbor_[0][0] << " " <<  neighbor_[0][1] << std::endl;
	//std::cout << "Voisin dir 1 : " << ipatch << " : " <<  neighbor_[1][0] << " " <<  neighbor_[1][1] << std::endl;
	

	min_local.resize(params.nDim_field, 0.);
	max_local.resize(params.nDim_field, 0.);
	cell_starting_global_index.resize(params.nDim_field, 0);
	for (int i = 0 ; i<params.nDim_field ; i++) {
	    min_local[i] = smpi->getDomainLocalMin(i) + Pcoordinates[i]*params.n_space[i]*params.cell_length[i];
	    max_local[i] = min_local[i] + params.n_space[i]*params.cell_length[i];
	    cell_starting_global_index[i] += Pcoordinates[i]*params.n_space[i];
	}

	std::cout << "Create patch\n\n";

	vecSpecies = SpeciesFactory::createVector(params, smpi, this);

	/* // + min_loc/cell_index(ref smpi,  & sort) // OK through this 
 * 	   std::cout << "Patch created\n";*/
	// + new n_space -> in PatchFactory
	// patchID : ok through coord
	// create Pos : OK

	// -> partBoundCond : min/max_loc (smpi)
	EMfields   = ElectroMagnFactory::create(params, laser_params, smpi, this);
	// + patchId + new n_space (now = params by smpi) + BC
	// -> Neighbors to define !!
	
	Interp     = InterpolatorFactory::create(params, smpi, this);               // + patchId -> idx_domain_begin (now = ref smpi)
	Proj       = ProjectorFactory::create(params, smpi, this);                  // + patchId -> idx_domain_begin (now = ref smpi)
	
};


//Functions for manipulating Hilbert curves in an arbitrary number of dimensions and arbitrary size.
//Taken from Chris Hamilton, Technical Report CS-2006-07, Faculty of computer Science, Halifax.

//Bitwise rotation operators.
//Works only for shift <= dim.
unsigned int Patch::rotl(unsigned int value, unsigned int shift, unsigned int dim) 
{
    //Evaluate left side of the rotated value then right side of it. Finally retain only the first dim bits.
    return ((value << shift) | (value >> (dim - shift))) & ((1<<dim) - 1);
}
unsigned int Patch::rotr(unsigned int value, unsigned int shift, unsigned int dim) 
{
    //Evaluate right side of the rotated value then left side of it. Finally retain only the first dim bits.
    return ((value >> shift) | (value << (dim - shift))) & ((1<<dim) - 1);
}

//Generates the binary reflected Gray Code
unsigned int Patch::gc(unsigned int i)
{
return i^(i>>1);
}
//Given a non-negative integer g, calculates the non-negative integer i such that gc(i)=g.
unsigned int Patch::gcinv(unsigned int g)
{
    unsigned int i,j;
    i=g;
    j=1;
    while ((1<<j) <= g){
    i = i ^ (g >> j);    
    j++;
    }
    return i;
}

// Tsb = trailing set bit. It is the number of trailing set bits in the binary representation of i.
// tsb is also the inter sub-hypercube directio, g(i).
unsigned int Patch::tsb(unsigned int i)
{
    unsigned int k;
    k = 0;
    while (i & 1) {
        i = i>>1;
        k++;
    }
    return k;
}
// Direction computes the sequence of intra sub-hypercube direction, d(i) for 0 <= i < 2^dim.
unsigned int Patch::direction(unsigned int i, unsigned int dim)
{
    if (i == 0) {
        return 0;
    }else if (i & 1){
        return tsb(i)%dim;
    }else {
        return tsb(i-1)%dim;
    }
}
// Entry computes the sequence of entry points, e(i) for 0 <= i < 2^dim.
    unsigned int Patch::entry(unsigned int i)
{
    if (i == 0) {
        return 0;
    }else {
        return gc(2*((i-1)/2));
    }
}
//!Ted is the transformation such that the gc ordering of sub-hypercubes in the Hilbert curve defined by e and d will map tot he standard binary reflected gc.
void Patch::ted(unsigned int e, unsigned int d, unsigned int *b, unsigned int dim)
{
    *b = rotr( *b ^ e, d+1, dim );
    return;
}
void Patch::tedinv(unsigned int e, unsigned int d, unsigned int *b, unsigned int dim)
{
    *b = rotl( *b , d+1, dim ) ^ e;
    return;
}
//!Hilbert index2D calculates the Hilbert index h of a patch of coordinates x,y for a simulation box with 2^m patches per side (2^(2*m) patches in total).
unsigned int Patch::hilbertindex(unsigned int m, unsigned int x, unsigned int y)
{
    unsigned int e,d,h,l,w;
    h = 0;
    e = 0;
    d = 0;
    for (int i = m-1; i>=0; i--){
        l = bit(y,i)*2 + bit(x,i); //ith bit of y at the leftmost position of l, and ith bit of x at the rightmost position of l.
        ted(e,d, &l, 2); 
        w = gcinv(l);
        e = e ^ (rotl(entry(w), d+1, 2));
        d = (d + direction(w, 2) +1 )%2 ;
        h = (h<<2)|w;
    }
    return h;
}
//!Hilbert index3D calculates the Hilbert index h of a patch of coordinates x,y,z for a simulation box with 2^m patches per side (2^(3*m) patches in total).
unsigned int Patch::hilbertindex(unsigned int m, unsigned int x, unsigned int y, unsigned int z)
{
    unsigned int e,d,h,l,w;
    h = 0;
    e = 0;
    d = 0;
    for (int i = m-1; i>=0; i--){
        l = bit(z,i)*4 + bit(y,i)*2 + bit(x,i); 
        ted(e,d, &l, 3); 
        w = gcinv(l);
        e = e ^ (rotl(entry(w), d+1, 3));
        d = (d + direction(w, 3) + 1)%3 ;
        h = (h<<3)|w;
    }
    return h;
}

//!Hilbert index2D inv  calculates the coordinates x,y of the patch of Hilbert index h in a simulation box with 2^m patches per side (2^(2*m) patches in total).
void Patch::hilbertindexinv(unsigned int m, unsigned int* x, unsigned int* y, unsigned int h)
{
    unsigned int e,d,l,w;
    e = 0;
    d = 0;
    *x=0;
    *y=0;
    for( int i = m-1; i>=0 ; i--){
        w = ((bit(h,2*i+1))<<1) + bit(h,2*i);
        l = gc(w);
        tedinv(e,d,&l,2);
        setbit(x,(unsigned int)i,bit(l,0));
        setbit(y,(unsigned int)i,bit(l,1));
        e = e ^ (rotl(entry(w), d+1, 2));
        d = (d + direction(w, 2) +1 )%2 ;
    }
    return;
}

//!extractMask extracts a mask µ indicating which axes are active at a given iteration i of the compact hilbert index.
//! For a simulation box with 2^m0 patches along X and 2^m1 patches along Y.
unsigned int Patch::extractmask(unsigned int m0,unsigned int  m1,unsigned int  i)
{
    unsigned int mu;
    mu = 0;
    if (m1 > i) mu = mu | 1;
    mu = mu << 1;
    if (m0 > i) mu = mu | 1;
    return mu;
}
//!extractMask extracts a mask µ indicating which axes are active at a given iteration i of the compact hilbert index.
//! For a simulation box with 2^m0 patches along X, 2^m1 patches along Y and 2^m2 patches along Z.
unsigned int Patch::extractmask(unsigned int m0,unsigned int  m1,unsigned int  m2,unsigned int  i)
{
    unsigned int mu;
    mu = 0;
    if (m2 > i) mu = mu | 1;
    mu = mu << 1;
    if (m1 > i) mu = mu | 1;
    mu = mu << 1;
    if (m0 > i) mu = mu | 1;
    return mu;
}
//!Gray Code Rank.
unsigned int Patch::gcr(unsigned int dim, unsigned int mu,unsigned int i)
{
    unsigned int r;
    r = 0;
    for (int k = dim-1; k>=0; k--){
        if( bit(mu, k) ) r = (r << 1) | bit(i,k);
    }
    return r;
}
//!Gray Code Rank Inverse.
unsigned int Patch::gcrinv(unsigned int dim, unsigned int mu,unsigned int pi, unsigned int r)
{
    unsigned int g,i,j;
    g = 0;
    i = 0;
    j = 0;
    for (unsigned int k=0; k < dim ; k++) j += bit(mu,k) ; //Counts the number of 1 in the binary representation of mu.
    j--; //At this point, j = ||mu|| - 1

    for (int k=dim-1; k >=0 ; k--){
        if(  bit(mu,k) ){
            setbit(&i, k, bit(r,j));
            setbit(&g, k, (bit(i,k) + bit(i,k+1))%2 );
            j--;
        } else {
            setbit(&g, k, bit(pi,k));
            setbit(&i, k, (bit(g,k)+bit(i,k+1))%2);
        }
    }
    return i;
}
//!Get kth bit of i.
unsigned int Patch::bit(unsigned int i, unsigned int k)
{
    return (i>>k)&1 ;
}
//!Set kth bit of i to value.
void Patch::setbit(unsigned int* i, unsigned int k, unsigned int value)
{
    *i = (*i & ~(1<<k)) | (value<<k);
    return;
}

//The "compact" version of the functions allows a different number of patch in each direction.

//!Compact Hilbert index2D calculates the compact Hilbert index h of a patch of coordinates x,y for a simulation box with 2^mi patches per side (2^(m0+m1)) patches in total).
unsigned int Patch::compacthilbertindex(unsigned int m0, unsigned int m1, unsigned int x, unsigned int y)
{
    unsigned int h,e,d,m,mu,w,l,r;
    h=0;
    e=0;
    d=0;
    m=std::max(m0, m1);
    for (int i = m-1; i>=0; i--){
        mu = extractmask(m0, m1, (unsigned int)i);
        mu = rotr(mu,d+1,2);
        l = bit(y,i)*2 + bit(x,i); 
        ted(e,d, &l, 2); 
        w = gcinv(l);
        r = gcr(2,mu,w);
        e = e ^ (rotl(entry(w), d+1, 2));
        d = (d + direction(w, 2) + 1)%2 ;
        for (unsigned int k=0; k < 2 ; k++) h = h << bit(mu,k) ;
        h = h | r ;
    }
    return h;
}
//!Compact Hilbert index3D calculates the compact Hilbert index h of a patch of coordinates x,y,z for a simulation box with 2^mi patches per side (2^(m0+m1+m2)) patches in total).
unsigned int Patch::compacthilbertindex(unsigned int m0, unsigned int m1, unsigned int m2, unsigned int x, unsigned int y, unsigned int z)
{
    unsigned int h,e,d,m,mu,w,l,r;
    h=0;
    e=0;
    d=0;
    m=std::max(std::max(m0, m1), m2);
    for (int i = m-1; i>=0; i--){
        mu = extractmask(m0, m1, m2, (unsigned int)i);
        mu = rotr(mu,d+1,3);
        l = bit(z,i)*4 + bit(y,i)*2 + bit(x,i); 
        ted(e,d, &l, 3); 
        w = gcinv(l);
        r = gcr(3,mu,w);
        e = e ^ (rotl(entry(w), d+1, 3));
        d = (d + direction(w, 3) + 1)%3 ;
        for (unsigned int k=0; k < 3 ; k++) h = h << bit(mu,k) ;
        h = h | r ;

    }
    return h;
}
//!Hilbert index inverse calculates the coordinates x,y of a patch for a given Hilbert index h in a simulation box with 2^mi patches per side (2^(m0+m1) patches in total)
//2D version
void Patch::compacthilbertindexinv(unsigned int m0, unsigned int m1, unsigned int* x, unsigned int* y, unsigned int h)
{
    unsigned int e,d,k,mu,l,r,mmax,msum,w,norm,pi;
    e=0;
    d=0;
    k=0;
    *x=0;
    *y=0;
    msum=m0+m1;
    mmax=std::max(m0, m1);
    for (int i=mmax-1; i>=0; i--){
        mu = extractmask(m0, m1, (unsigned int)i);
        mu = rotr(mu,d+1,2);
        pi = rotr(e,d+1,2) & ~mu;
        norm=0;
        for (unsigned int n=0; n < 2 ; n++) norm += bit(mu,n) ;
        r=0;
        for (unsigned int n=0; n < norm ; n++){
            r = r << 1;
            r += bit(h,msum-1-k-n) ;
        }
        k += norm;
        w = gcrinv(2,mu,pi,r);
        l = gc(w);
        tedinv(e,d,&l,2);
        setbit(x,(unsigned int)i,bit(l,0));
        setbit(y,(unsigned int)i,bit(l,1));
        e = e ^ (rotl(entry(w), d+1, 2));
        d = (d + direction(w, 2) + 1)%2 ;
    } 
    return;
}
//3D version
//Coordinates x,y,z of a patch of index h in a simulation box with a total number of patch = 2^(m0+m1+m2)
void Patch::compacthilbertindexinv(unsigned int m0, unsigned int m1, unsigned int m2,  unsigned int* x, unsigned int* y, unsigned int* z, unsigned int h)
{
    unsigned int e,d,k,mu,l,r,mmax,msum,w,norm,pi;
    e=0;
    d=0;
    k=0;
    *x=0;
    *y=0;
    *z=0;
    msum=m0+m1+m2;
    mmax=std::max(std::max(m0, m1), m2);
    for (int i=mmax-1; i>=0; i--){
        mu = extractmask(m0, m1, m2, (unsigned int)i);
        mu = rotr(mu,d+1,3);
        pi = rotr(e,d+1,3) & ~mu;
        norm=0;
        for (unsigned int n=0; n < 3 ; n++) norm += bit(mu,n) ;
        r=0;
        for (unsigned int n=0; n < norm ; n++){
            r = r << 1;
            r += bit(h,msum-1-k-n) ;
        }
        k += norm;
        w = gcrinv(3,mu,pi,r);
        l = gc(w);
        tedinv(e,d,&l,3);
        setbit(x,(unsigned int)i,bit(l,0));
        setbit(y,(unsigned int)i,bit(l,1));
        setbit(z,(unsigned int)i,bit(l,2));
        e = e ^ (rotl(entry(w), d+1, 3));
        d = (d + direction(w, 3) + 1)%3 ;
    } 
    return;
}

void Patch::dynamics(double time_dual, SmileiMPI *smpi, PicParams &params, SimWindow* simWindow, int diag_flag)
{
    for (unsigned int ispec=0 ; ispec<params.n_species; ispec++) {
	if ( vecSpecies[ispec]->isProj(time_dual, simWindow) || diag_flag  ){    
	    vecSpecies[ispec]->dynamics(time_dual, ispec, EMfields, Interp, Proj, smpi, params, simWindow, diag_flag);
	}
    }

}


void Patch::exchParticles(SmileiMPI* smpi, int ispec, PicParams &params, int tid, int iDim)
{
    Particles &cuParticles = (*vecSpecies[ispec]->particles);
    std::vector<int>* cubmin = &vecSpecies[ispec]->bmin;
    std::vector<int>* cubmax = &vecSpecies[ispec]->bmax;
    
    std::vector< std::vector<int> >* indexes_of_particles_to_exchange_per_thd = &vecSpecies[ispec]->indexes_of_particles_to_exchange_per_thd;
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


        Particles diagonalParticles;
        diagonalParticles.initialize(0,cuParticles.dimension());     

	// A définir : buff_index_send
	int nbNeighbors_(2);
	std::vector<int> buff_index_send[nbNeighbors_];
	int buff_index_recv_sz[nbNeighbors_];
	for (int i=0 ; i<nbNeighbors_ ; i++) {
	    buff_index_send[i].resize(0);
	    buff_index_recv_sz[i] = 0;
	}    
	// A définir : buff_index_send

        for (int i=0 ; i<n_part_send ; i++) {
            iPart = indexes_of_particles_to_exchange[i];
	    if ( cuParticles.position(iDim,iPart) < min_local[iDim] ) {
		if (neighbor_[iDim][0]!=-2)
		    buff_index_send[0].push_back( indexes_of_particles_to_exchange[i] );
		else // To another MPI process
		    cuParticles.cp_particle(indexes_of_particles_to_exchange[i], smpi->interParticles);
	    }
	    else if ( ( cuParticles.position(iDim,iPart) >= max_local[iDim]) ) {
		if (neighbor_[iDim][1]!=-2)
		    buff_index_send[1].push_back( indexes_of_particles_to_exchange[i] );
	        else // To another MPI process
		    cuParticles.cp_particle(indexes_of_particles_to_exchange[i], smpi->interParticles);
	    }
	    else if ( !(cuParticles.is_part_in_domain(iPart, this) ) ) {
		// at the end of exchangeParticles, diagonalParticles will be reinjected 
		// at the end of cuParticles & indexes_of_particles_to_exchange_per_thd[0] for next iDim
		cuParticles.cp_particle(indexes_of_particles_to_exchange[i], diagonalParticles);
	    }
	    else { // particle will be deleted, if supp_particle particles still in the domain
	    }
        } // END for iPart = f(i)

        
        Particles partVectorSend[2];
        partVectorSend[0].initialize(0,cuParticles.dimension());
        partVectorSend[1].initialize(0,cuParticles.dimension());
        Particles partVectorRecv[2];
        partVectorRecv[0].initialize(0,cuParticles.dimension());
        partVectorRecv[1].initialize(0,cuParticles.dimension());
        
        /********************************************************************************/
        // Exchange particles
        /********************************************************************************/
            
	MPI_Status sstat    [2];
	MPI_Status rstat    [2];
	MPI_Request srequest[2];
	MPI_Request rrequest[2];
            

	/********************************************************************************/
	// Exchange number of particles to exchange to establish or not a communication
	/********************************************************************************/
	for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
	    if (neighbor_[iDim][iNeighbor]!=-2) {
		n_part_send = (buff_index_send[iNeighbor]).size();
		// tag = sign//hindex//neighbor_[iDim][iNeighbor]
		stringstream stag("");
		stag << hindex << "00" << neighbor_[iDim][iNeighbor];
		int tag(0);
		stag >> tag;
		MPI_Isend( &n_part_send, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &(srequest[iNeighbor]) );
	    } // END of Send
	    else
	        n_part_send = 0;
	    if (neighbor_[iDim][(iNeighbor+1)%2]!=-2) {
		buff_index_recv_sz[(iNeighbor+1)%2] = 0;
		stringstream stag("");
		stag << neighbor_[iDim][iNeighbor] << "00" << hindex;
		int tag(0);
		stag >> tag;
		MPI_Irecv( &(buff_index_recv_sz[(iNeighbor+1)%2]), 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &(rrequest[(iNeighbor+1)%2]) );
	    }
	    else 
	        buff_index_recv_sz[(iNeighbor+1)%2] = 0;
	}
	//smpi->barrier();
            
	/********************************************************************************/
	// Wait for end of communications over number of particles
	/********************************************************************************/
	for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
	    if (neighbor_[iDim][iNeighbor]!=-2) {
		MPI_Wait( &(srequest[iNeighbor]), &(sstat[iNeighbor]) );
	    }
	    if (neighbor_[iDim][(iNeighbor+1)%2]!=-2) {
		MPI_Wait( &(rrequest[(iNeighbor+1)%2]), &(rstat[(iNeighbor+1)%2]) );
		if (buff_index_recv_sz[(iNeighbor+1)%2]!=0) {
		    partVectorRecv[(iNeighbor+1)%2].initialize( buff_index_recv_sz[(iNeighbor+1)%2], cuParticles.dimension());
		}
	    }
	}
	//smpi->barrier();

#ifdef _PATATE
            
	/********************************************************************************/
	// Define buffers to exchange buff_index_send[iNeighbor].size();
	/********************************************************************************/
            
            
	/********************************************************************************/
	// Proceed to effective Particles' communications
	/********************************************************************************/

	// Number of properties per particles = nDim_Particles + 3 + 1 + 1
	int nbrOfProp( 7 );
	MPI_Datatype typePartSend, typePartRecv;

	for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
                
	    // n_part_send : number of particles to send to current neighbor
	    n_part_send = (buff_index_send[iNeighbor]).size();
	    if ( (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) && (n_part_send!=0) ) {
		double x_max = params.cell_length[iDim]*( params.n_space_global[iDim] );
		for (int iPart=0 ; iPart<n_part_send ; iPart++) {
		    if (periods_[iDim]==1) {
			// Enabled periodicity
			if ( ( iNeighbor==0 ) &&  (coords_[iDim] == 0 ) &&( cuParticles.position(iDim,buff_index_send[iNeighbor][iPart]) < 0. ) ) {
			    cuParticles.position(iDim,buff_index_send[iNeighbor][iPart])     += x_max;
			}
			else if ( ( iNeighbor==1 ) &&  (coords_[iDim] == number_of_procs[iDim]-1 ) && ( cuParticles.position(iDim,buff_index_send[iNeighbor][iPart]) >= x_max ) ) {
			    cuParticles.position(iDim,buff_index_send[iNeighbor][iPart])     -= x_max;
			}
		    }
		    cuParticles.cp_particle(buff_index_send[iNeighbor][iPart], partVectorSend[iNeighbor]);
		}

		typePartSend = createMPIparticles( &(partVectorSend[iNeighbor]), nbrOfProp );
		MPI_Isend( &((partVectorSend[iNeighbor]).position(0,0)), 1, typePartSend, neighbor_[iDim][iNeighbor], 0, SMILEI_COMM_2D, &(srequest[iNeighbor]) );
		MPI_Type_free( &typePartSend );

	    } // END of Send
                
	    n_part_recv = buff_index_recv_sz[(iNeighbor+1)%2];
	    if ( (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
		typePartRecv = createMPIparticles( &(partVectorRecv[(iNeighbor+1)%2]), nbrOfProp );
		MPI_Irecv( &((partVectorRecv[(iNeighbor+1)%2]).position(0,0)), 1, typePartRecv,  neighbor_[iDim][(iNeighbor+1)%2], 0, SMILEI_COMM_2D, &(rrequest[(iNeighbor+1)%2]) );
		MPI_Type_free( &typePartRecv );

	    } // END of Recv
                
	} // END for iNeighbor
            
            
	/********************************************************************************/
	// Wait for end of communications over Particles
	/********************************************************************************/
	for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
                
	    n_part_send = buff_index_send[iNeighbor].size();
	    n_part_recv = buff_index_recv_sz[(iNeighbor+1)%2];
                
	    if ( (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) && (n_part_send!=0) ) {
		MPI_Wait( &(srequest[iNeighbor]), &(sstat[iNeighbor]) );
	    }
                
	    if ( (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
		MPI_Wait( &(rrequest[(iNeighbor+1)%2]), &(rstat[(iNeighbor+1)%2]) );
                    
		// Extract corner particles, not managed in the following process 
		// but reinjected at the end in the main list
		for (int iPart=0 ; iPart<n_part_recv; iPart++ )
		    if ( !(partVectorRecv[(iNeighbor+1)%2]).is_part_in_domain(iPart, this) )
			(partVectorRecv[(iNeighbor+1)%2]).cp_particle(iPart, diagonalParticles);
		for (int iPart=n_part_recv-1 ; iPart>=0; iPart-- ) {
		    if ( !(partVectorRecv[(iNeighbor+1)%2]).is_part_in_domain(iPart, this) ) {
			(partVectorRecv[(iNeighbor+1)%2]).erase_particle(iPart);
			buff_index_recv_sz[(iNeighbor+1)%2]--;
		    }
		}
                    
	    }
                
	}
	barrier();
	/********************************************************************************/
	// Clean lists of indexes of particle to exchange per neighbor
	/********************************************************************************/
	for (int i=0 ; i<nbNeighbors_ ; i++)
	    buff_index_send[i].clear();
            
        
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
	if (iDim==1) {
	    //1) Count particles coming from south and north
	    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
		n_part_recv = buff_index_recv_sz[iNeighbor];
		for (unsigned int j=0; j<n_part_recv ;j++){
		    ii = int((partVectorRecv[iNeighbor].position(0,j)-min_local[0])/dbin);//bin in which the particle goes.
		    shift[ii+1]++; // It makes the next bins shift.
		}
	    }
	}
	if (iDim==0) {
	    //2) Add particles coming from west and east
	    shift[1] += buff_index_recv_sz[0];//Particles coming from south all go to bin 0 and shift all the other bins.
	    shift[(*cubmax).size()] += buff_index_recv_sz[1];//Used only to count the total number of particles arrived.
	}
        
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
	if (iDim==0) {
	    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
		n_part_recv = buff_index_recv_sz[iNeighbor];
		if ( (neighbor_[0][iNeighbor]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
		    ii = iNeighbor*((*cubmax).size()-1);//0 if iNeighbor=0(particles coming from West) and (*cubmax).size()-1 otherwise.
		    partVectorRecv[iNeighbor].overwrite_part2D(0, cuParticles,(*cubmax)[ii],n_part_recv);
		    (*cubmax)[ii] += n_part_recv ;
		}
	    }
	}
        //iDim == 1) is the difficult case, when particles can arrive in any bin.
	if (iDim==1) {
	    for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
		n_part_recv = buff_index_recv_sz[iNeighbor];
		if ( (neighbor_[1][iNeighbor]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
		    for(unsigned int j=0; j<n_part_recv; j++){
			ii = int((partVectorRecv[iNeighbor].position(0,j)-min_local[0])/dbin);//bin in which the particle goes.
			partVectorRecv[iNeighbor].overwrite_part2D(j, cuParticles,(*cubmax)[ii]);
			(*cubmax)[ii] ++ ;
		    }
		}
	    }
	}
#endif
        
        // Inject corner particles at the end of the list, update bmax
	//if (iDim==cuParticles.dimension()-1) cout << "Number of diag particles " << diagonalParticles.size() << endl;
        for (int iPart = 0 ; iPart<diagonalParticles.size() ; iPart++) {
            diagonalParticles.cp_particle(iPart, cuParticles);
            (*indexes_of_particles_to_exchange_per_thd)[0].push_back(cuParticles.size()-1);
            (*cubmax)[(*cubmax).size()-1]++;
        }
        
    }//end of omp master

}
