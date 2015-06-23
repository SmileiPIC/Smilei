#include "Patch.h"
#include <cstdlib>
#include <iostream>
#include <iomanip>

using namespace std;

//int buildtag(int send, int recv);

Patch::Patch(PicParams& params, DiagParams &diag_params, LaserParams& laser_params, SmileiMPI* smpi, unsigned int m0, unsigned int m1, unsigned int m2, unsigned int ipatch) {


//Neighborhood definition in 2D:
//   
//   Y axis
//   ^
//   |   6       7          8
//   |   3    4(self)       5
//       0       1          2    --> X axis

        int xcall, ycall;
        hindex = ipatch;
	std::cout << smpi->getRank() << ", mypatch is : " << hindex << std::endl;        
        
        if ( params.geometry == "1d3v" ) {
            mi.resize(1);
            Pcoordinates.resize(1);
            mi[0] = m0;
            Pcoordinates[0] = hindex;
	    MPI_neighborhood_.resize(3);
	    patch_neighborhood_.resize(3);
        }
        else if ( params.geometry == "2d3v" ) {
            mi.resize(2);
            Pcoordinates.resize(2);
            mi[0] = m0;
            mi[1] = m1;
            generalhilbertindexinv(m0, m1, &Pcoordinates[0], &Pcoordinates[1], hindex);

	    /////////////////////////////////////
	    //  Define local domain
	    /////////////////////////////////////

	    MPI_neighborhood_.resize(9);
	    patch_neighborhood_.resize(9);
        }
        else {
            mi.resize(3);
            Pcoordinates.resize(3);
            mi[0] = m0;
            mi[1] = m1;
            mi[2] = m2;
            generalhilbertindexinv(m0, m1, m2, &Pcoordinates[0], &Pcoordinates[1], &Pcoordinates[2], hindex);
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
	if (params.bc_em_type_long=="periodic") xcall = xcall%((1<<m0));
	neighbor_[0][0] = generalhilbertindex( m0, m1, xcall, ycall);
#ifdef _PATCH_DEBUG
        cout << xcall << " " << ycall << " " << neighbor_[0][0] << endl;
#endif
        xcall = Pcoordinates[0]+1;
	if (params.bc_em_type_long=="periodic") xcall = xcall%((1<<m0));
	neighbor_[0][1] = generalhilbertindex( m0, m1, xcall, ycall);
        xcall = Pcoordinates[0];
        ycall = Pcoordinates[1]-1;
	if (params.bc_em_type_trans=="periodic") ycall = ycall%((1<<m1));
	neighbor_[1][0] = generalhilbertindex( m0, m1, xcall, ycall);
        ycall = Pcoordinates[1]+1;
	if (params.bc_em_type_trans=="periodic") ycall = ycall%((1<<m1));
	neighbor_[1][1] = generalhilbertindex( m0, m1, xcall, ycall);


        xcall = Pcoordinates[0]+1;
	if (params.bc_em_type_long=="periodic") xcall = xcall%((1<<m0));
	corner_neighbor_[1][1] = generalhilbertindex( m0, m1, xcall, ycall);
        xcall = Pcoordinates[0]-1;
	if (params.bc_em_type_long=="periodic") xcall = xcall%((1<<m0));
	corner_neighbor_[0][1] = generalhilbertindex( m0, m1, xcall, ycall);
        ycall = Pcoordinates[1]-1;
	if (params.bc_em_type_trans=="periodic") ycall = ycall%((1<<m1));
	corner_neighbor_[0][0] = generalhilbertindex( m0, m1, xcall, ycall);
        xcall = Pcoordinates[0]+1;
	if (params.bc_em_type_long=="periodic") xcall = xcall%((1<<m0));
	corner_neighbor_[1][0] = generalhilbertindex( m0, m1, xcall, ycall);


        patch_neighborhood_[0] = corner_neighbor_[0][0];
        patch_neighborhood_[1] = neighbor_[1][0];
        patch_neighborhood_[2] = corner_neighbor_[1][0];
        patch_neighborhood_[3] = neighbor_[0][0];
        patch_neighborhood_[4] = hindex;
        patch_neighborhood_[5] = neighbor_[0][1];
        patch_neighborhood_[6] = corner_neighbor_[0][1];
        patch_neighborhood_[7] = neighbor_[1][1];
        patch_neighborhood_[8] = corner_neighbor_[1][1];



	for ( int z = 0 ; z < 1+2*(params.nDim_field == 3) ; z++ ) {
	    for ( int y = 0 ; y < 1+2*(params.nDim_field >= 2) ; y++ ) {
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

	vecSpecies = SpeciesFactory::createVector(params, smpi, this);

	/* // + min_loc/cell_index(ref smpi,  & sort) // OK through this */
	// + new n_space -> in PatchFactory
	// patchID : ok through coord
	// create Pos : OK

	// -> partBoundCond : min/max_loc (smpi)
	EMfields   = ElectroMagnFactory::create(params, laser_params, smpi, this);
	// + patchId + new n_space (now = params by smpi) + BC
	// -> Neighbors to define !!
	
	Interp     = InterpolatorFactory::create(params, smpi, this);               // + patchId -> idx_domain_begin (now = ref smpi)
	Proj       = ProjectorFactory::create(params, smpi, this);                  // + patchId -> idx_domain_begin (now = ref smpi)

	Diags = new Diagnostic(params,diag_params, smpi, this);
	
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
unsigned int Patch::hilbertindex(unsigned int m, unsigned int x, unsigned int y, unsigned int *einit, unsigned int *dinit)
{
    unsigned int e,d,h,l,w;
    e = *einit;
    d = *dinit;
    h = 0;
    for (int i = m-1; i>=0; i--){
        l = bit(y,i)*2 + bit(x,i); //ith bit of y at the leftmost position of l, and ith bit of x at the rightmost position of l.
        ted(e,d, &l, 2); 
        w = gcinv(l);
        e = e ^ (rotl(entry(w), d+1, 2));
        d = (d + direction(w, 2) +1 )%2 ;
        h = (h<<2)|w;
    }

    *einit = e;
    *dinit = d;

    return h;
}
//!Hilbert index3D calculates the Hilbert index h of a patch of coordinates x,y,z for a simulation box with 2^m patches per side (2^(3*m) patches in total).
unsigned int Patch::hilbertindex(unsigned int m, unsigned int x, unsigned int y, unsigned int z,unsigned int einit, unsigned int dinit)
{
    unsigned int e,d,h,l,w;
    h = 0;
    e = einit;
    d = dinit;
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
void Patch::hilbertindexinv(unsigned int m, unsigned int* x, unsigned int* y, unsigned int h, unsigned int einit, unsigned int dinit)
{
    unsigned int e,d,l,w;
    e = einit;
    d = dinit;
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
void Patch::hilbertindexinv(unsigned int m, unsigned int* x, unsigned int* y, unsigned int* z, unsigned int h, unsigned int einit, unsigned int dinit)
{
    unsigned int e,d,l,w;
    e = einit;
    d = dinit;
    *x=0;
    *y=0;
    *z=0;
    for( int i = m-1; i>=0 ; i--){
        w = ((bit(h,3*i+2))<<2) + ((bit(h,3*i+1))<<1) + bit(h,3*i);
        l = gc(w);
        tedinv(e,d,&l,3);
        setbit(x,(unsigned int)i,bit(l,0));
        setbit(y,(unsigned int)i,bit(l,1));
        setbit(z,(unsigned int)i,bit(l,2));
        e = e ^ (rotl(entry(w), d+1, 3));
        d = (d + direction(w, 3) +1 )%3 ;
    }
    return;
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

//The "general" versions of the functions allow a different number of patch in each direction.

//!General Hilbert index2D calculates the  Hilbert index h of a patch of coordinates x,y for a simulation box with 2^mi patches per side (2^(m0+m1)) patches in total).
unsigned int Patch::generalhilbertindex(unsigned int m0, unsigned int m1, unsigned int x, unsigned int y, unsigned int *einit, unsigned int *dinit)
{
    if(x%((1<<m0)) != x || y%((1<<m1)) != y ) return MPI_PROC_NULL ;

    unsigned int h,mmin,mmax,l,localx,localy,*target;
    h=0;
    *dinit=0;
    localx = x;
    localy = y;
    if (m0 >= m1){
       target = &localx;
       mmin = m1;
       mmax = m0;
    } else {
       target = &localy;
       *dinit=1;
       mmin = m0;
       mmax = m1;
    }
    for (int i= mmax-1; i >= mmin ; i--){
       l = bit(*target,i); 
       h += l*(1<<(i+mmin));
       *target -= l*(1<<i);
    }
    if (mmin > 0) {
        h += hilbertindex(mmin,localx,localy,einit,dinit);
    }
return h;
}
unsigned int Patch::generalhilbertindex(unsigned int m0, unsigned int m1, unsigned int x, unsigned int y)
{
    if(x%((1<<m0)) != x || y%((1<<m1)) != y ) return MPI_PROC_NULL ;

    unsigned int h,mmin,mmax,l,localx,localy,*target,einit,dinit;
    h=0;
    dinit=0;
    einit=0;
    localx = x;
    localy = y;
    if (m0 >= m1){
       target = &localx;
       mmin = m1;
       mmax = m0;
    } else {
       target = &localy;
       dinit=1;
       mmin = m0;
       mmax = m1;
    }
    for (int i= mmax-1; i >= mmin ; i--){
       l = bit(*target,i); 
       h += l*(1<<(i+mmin));
       *target -= l*(1<<i);
    }
    if (mmin > 0) {
        h += hilbertindex(mmin,localx,localy,&einit,&dinit);
    }
return h;
}

//!General Hilbert index3D calculates the compact Hilbert index h of a patch of coordinates x,y,z for a simulation box with 2^mi patches per side (2^(m0+m1+m2)) patches in total).
unsigned int Patch::generalhilbertindex(unsigned int m0, unsigned int m1, unsigned int m2, unsigned int x, unsigned int y, unsigned int z)
{
    if(x%((1<<m0)) != x || y%((1<<m1)) != y || z%((1<<m2)) != z) return MPI_PROC_NULL ;

    unsigned int h,e,d,*einit,*dinit,dimmin,dimmax,dimmed,l,localx,localy,localz, mi[3],localp[3],tempp[3],mmin;
    h=0;
    e=0;
    d=0;
    dinit=&d;
    einit=&e;
    //Store positions and dimensions in arrays
    localp[0] = x;
    localp[1] = y;
    localp[2] = z;
    mi[0] = m0;
    mi[1] = m1;
    mi[2] = m2;
    //Compare dimension sizes
    if ((m0 >= m1) && (m0 >= m2) ){
        dimmax = 0;
    } else if ((m1 > m0) && (m1 >= m2)){
        dimmax = 1;
    } else {
        dimmax = 2;
    }
    if (mi[(dimmax+1)%3] >= mi[(dimmax+2)%3]){
        dimmed = (dimmax+1)%3;
        dimmin = (dimmax+2)%3;
    } else {
        dimmed = (dimmax+2)%3;
        dimmin = (dimmax+1)%3;
    }
    // First approach on a flattened 2D grid along dimmax and dimmed. The 3D grid is projected along dimmin axis.
    tempp[dimmax] = localp[dimmax] >> mi[dimmin]; //Erase last mi[dimmin] bits. Not relevent for this phase.
    tempp[dimmed] = localp[dimmed] >> mi[dimmin]; //Erase last mi[dimmin] bits. Not relevent for this phase.
    h += generalhilbertindex(mi[dimmax]-mi[dimmin],mi[dimmed]-mi[dimmin],tempp[dimmax],tempp[dimmed],einit,dinit)*(1<<(3*mi[dimmin]));

    //Now in a local cube of side mi[dimmin]. The local entry point "einit" and initial direction "dinit" of the local hilbert curve has been
    //determined by the previous call to compacthilbertindex2.
    //Relative position in the local cube is given by the last mi[dimmin] bits of the position.
    tempp[dimmax] = localp[dimmax] & ((1<<mi[dimmin])-1); //Only keep the last mi[dimmin] bits.
    tempp[dimmed] = localp[dimmed] & ((1<<mi[dimmin])-1); //Only keep the last mi[dimmin] bits.
    tempp[dimmin] = localp[dimmin] & ((1<<mi[dimmin])-1); //Only keep the last mi[dimmin] bits.
    //Add local index to the previously calculated one.
    h += hilbertindex(mi[dimmin],tempp[dimmax], tempp[dimmed], tempp[dimmin], *einit, *dinit);
    


return h;
}
//!General Hilbert index inverse calculates the coordinates x,y of a patch for a given Hilbert index h in a simulation box with 2^mi patches per side (2^(m0+m1) patches in total)
//2D version
void Patch::generalhilbertindexinv(unsigned int m0, unsigned int m1, unsigned int* x, unsigned int* y, unsigned int h)
{
    unsigned int einit, dinit,l,localh, *target, shift ;
    int mmin,mmax;
    einit = 0;
    dinit = 0;
    shift = 0;
    localh = h;
    //Compare dimensions. Target points at the dimension which must be shifted. 
    if (m0 >= m1){
       target = x;
       mmin = m1;
       mmax = m0;
    } else {
       target = y;
       dinit=1;
       mmin = m0;
       mmax = m1;
    }
    //First define in which sub-hypercube of side 2^mmin the point is.
    for (int i= mmax+mmin-1; i >= mmin+mmin ; i--){
       l = bit(localh,i); 
       shift += l*(1<<(i-mmin));
       localh -= l*(1<<i);
    }
    //Run the cubic inversion algorithm in the sub hypercube.
    hilbertindexinv(mmin, x, y, localh, einit, dinit);
    //Shift the appropriate coordinate by the necessary value.
    *target += shift;

    return;

}
//3D version
void Patch::generalhilbertindexinv(unsigned int m0, unsigned int m1, unsigned int m2,  unsigned int* x, unsigned int* y, unsigned int* z, unsigned int h)
{
    unsigned int e,d,dimmin,dimmax,dimmed,l,localx,localy,localz, mi[3],*localp[3],tempp[3],localh;
    e=0;
    d=0;
    //Store positions and dimensions in arrays
    localp[0] = x;
    localp[1] = y;
    localp[2] = z;
    mi[0] = m0;
    mi[1] = m1;
    mi[2] = m2;
    //Compare dimension sizes
    if ((m0 >= m1) && (m0 >= m2) ){
        dimmax = 0;
    } else if ((m1 > m0) && (m1 >= m2)){
        dimmax = 1;
    } else {
        dimmax = 2;
    }
    if (mi[(dimmax+1)%3] >= mi[(dimmax+2)%3]){
        dimmed = (dimmax+1)%3;
        dimmin = (dimmax+2)%3;
    } else {
        dimmed = (dimmax+2)%3;
        dimmin = (dimmax+1)%3;
    }

    
    //Localize in which sub hypercube the point is. Do not account for the first 3*dimmin bits of h.
    localh = (h >> (mi[dimmin]*3));
    //Run the 2D inversion algorithm on the reduced domain.
    generalhilbertindexinv(mi[dimmax]-mi[dimmin],mi[dimmed]-mi[dimmin],localp[dimmax],localp[dimmed],localh);
    // Now local P stores the position of the cube in the 2D domain
    // We need to run the 3D inversion algorithm on this cube with the correct entry point and direction.
    // Run the 2D indexgenerator in order to evaluate e and d.
    localh = generalhilbertindex(mi[dimmax]-mi[dimmin],mi[dimmed]-mi[dimmin],*localp[dimmax],*localp[dimmed],&e,&d);
    //Transform coordinates in the global frame.
    *localp[dimmax] *= (1<<mi[dimmin]);
    *localp[dimmed] *= (1<<mi[dimmin]);
    *localp[dimmin] = 0 ;
    //Use only first bits of h for the local hypercube.
    localh = h & ((1<<(mi[dimmin]*3))-1);
    //Run the cubic inversion algorithm in the local sub hypercube.
    hilbertindexinv(mi[dimmin], &tempp[dimmax], &tempp[dimmed], &tempp[dimmin], localh, e, d);
    //Add results to the coordinates.
    *localp[dimmax] += tempp[dimmax];
    *localp[dimmed] += tempp[dimmed];
    *localp[dimmin] += tempp[dimmin];

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
}

//Old implmeentation inspired by Chris Hamilton. Not used anymore because the generated Hilbert curves are not continuous.
//void Patch::compacthilbertindexinv(unsigned int m0, unsigned int m1, unsigned int m2,  unsigned int* x, unsigned int* y, unsigned int* z, unsigned int h)
//{
//    unsigned int e,d,k,mu,l,r,mmax,msum,w,norm,pi;
//    e=0;
//    d=0;
//    k=0;
//    *x=0;
//    *y=0;
//    *z=0;
//    msum=m0+m1+m2;
//    mmax=std::max(std::max(m0, m1), m2);
//    for (int i=mmax-1; i>=0; i--){
//        mu = extractmask(m0, m1, m2, i);
//        mu = rotr(mu,d+1,3);
//        pi = rotr(e,d+1,3) & ~mu;
//        norm=0;
//        for (unsigned int n=0; n < 3 ; n++) norm += bit(mu,n) ;
//        r=0;
//        for (unsigned int n=0; n < norm ; n++){
//            r = r << 1;
//            r += bit(h,msum-1-k-n) ;
//        }
//        k += norm;
//        w = gcrinv(3,mu,pi,r);
//        l = gc(w);
//        tedinv(e,d,&l,3);
//        setbit(x,(unsigned int)i,bit(l,0));
//        setbit(y,(unsigned int)i,bit(l,1));
//        setbit(z,(unsigned int)i,bit(l,2));
//        e = e ^ (rotl(entry(w), d+1, 3));
//        d = (d + direction(w, 3) + 1)%3 ;
//    } 
//    return;
//}
//void Patch::compacthilbertindexinv(unsigned int m0, unsigned int m1, unsigned int* x, unsigned int* y, unsigned int h)
//{
//    unsigned int e,d,k,mu,l,r,mmax,msum,w,norm,pi;
//    e=0;
//    d=0;
//    k=0;
//    *x=0;
//    *y=0;
//    msum=m0+m1;
//    mmax=std::max(m0, m1);
//    for (int i=mmax-1; i>=0; i--){
//        mu = extractmask(m0, m1, i);
//        mu = rotr(mu,d+1,2);
//        pi = rotr(e,d+1,2) & ~mu;
//        norm=0;
//        for (unsigned int n=0; n < 2 ; n++) norm += bit(mu,n) ;
//        r=0;
//        for (unsigned int n=0; n < norm ; n++){
//            r = r << 1;
//            r += bit(h,msum-1-k-n) ;
//        }
//        k += norm;
//        w = gcrinv(2,mu,pi,r);
//        l = gc(w);
//        tedinv(e,d,&l,2);
//        setbit(x,(unsigned int)i,bit(l,0));
//        setbit(y,(unsigned int)i,bit(l,1));
//        e = e ^ (rotl(entry(w), d+1, 2));
//        d = (d + direction(w, 2) + 1)%2 ;
//    } 
//    return;
//}
//unsigned int Patch::compacthilbertindex(unsigned int m0, unsigned int m1, unsigned int m2, unsigned int x, unsigned int y, unsigned int z)
//{
//    unsigned int h,e,d,m,mu,w,l,r;
//    h=0;
//    e=0;
//    d=0;
//    m=std::max(std::max(m0, m1), m2);
//    for (int i = m-1; i>=0; i--){
//        mu = extractmask(m0, m1, m2, i);
//        mu = rotr(mu,d+1,3);
//        l = bit(z,i)*4 + bit(y,i)*2 + bit(x,i); 
//        ted(e,d, &l, 3); 
//        w = gcinv(l);
//        r = gcr(3,mu,w);
//        e = e ^ (rotl(entry(w), d+1, 3));
//        d = (d + direction(w, 3) + 1)%3 ;
//        for (unsigned int k=0; k < 3 ; k++) h = h << bit(mu,k) ;
//        h = h | r ;
//
//    }
//    return h;
//}
//unsigned int Patch::compacthilbertindex(unsigned int m0, unsigned int m1, unsigned int x, unsigned int y)
//{
//    unsigned int h,e,d,m,mu,w,l,r;
//    h=0;
//    e=0;
//    d=0;
//    m=std::max(m0, m1);
//    for (int i = m-1; i>=0; i--){
//        mu = extractmask(m0, m1, i);
//        mu = rotr(mu,d+1,2);
//        l = bit(y,i)*2 + bit(x,i); 
//        ted(e,d, &l, 2); 
//        w = gcinv(l);
//        r = gcr(2,mu,w);
//        e = e ^ (rotl(entry(w), d+1, 2));
//        d = (d + direction(w, 2) + 1)%2 ;
//        for (unsigned int k=0; k < 2 ; k++) h = h << bit(mu,k) ;
//        h = h | r ;
//    }
//    return h;
//}
//!Gray Code Rank Inverse.
//unsigned int Patch::gcrinv(unsigned int dim, unsigned int mu,unsigned int pi, unsigned int r)
//{
//    unsigned int g,i,j;
//    g = 0;
//    i = 0;
//    j = 0;
//    for (unsigned int k=0; k < dim ; k++) j += bit(mu,k) ; //Counts the number of 1 in the binary representation of mu.
//    j--; //At this point, j = ||mu|| - 1
//
//    for (int k=dim-1; k >=0 ; k--){
//        if(  bit(mu,k) ){
//            setbit(&i, k, bit(r,j));
//            setbit(&g, k, (bit(i,k) + bit(i,k+1))%2 );
//            j--;
//        } else {
//            setbit(&g, k, bit(pi,k));
//            setbit(&i, k, (bit(g,k)+bit(i,k+1))%2);
//        }
//    }
//    return i;
//}
//!extractMask extracts a mask µ indicating which axes are active at a given iteration i of the compact hilbert index.
//! For a simulation box with 2^m0 patches along X and 2^m1 patches along Y.
//unsigned int Patch::extractmask(unsigned int m0,unsigned int  m1, int  i)
//{
//    unsigned int mu;
//    mu = 0;
//    if (m1 > i) mu = mu | 1;
//    mu = mu << 1;
//    if (m0 > i) mu = mu | 1;
//    return mu;
//}
//!extractMask extracts a mask µ indicating which axes are active at a given iteration i of the compact hilbert index.
//! For a simulation box with 2^m0 patches along X, 2^m1 patches along Y and 2^m2 patches along Z.
//unsigned int Patch::extractmask(unsigned int m0,unsigned int  m1,unsigned int  m2, int  i)
//{
//    unsigned int mu;
//    mu = 0;
//    if (m2 > i) mu = mu | 1;
//    mu = mu << 1;
//    if (m1 > i) mu = mu | 1;
//    mu = mu << 1;
//    if (m0 > i) mu = mu | 1;
//    return mu;
//}
//!Gray Code Rank.
//unsigned int Patch::gcr(unsigned int dim, unsigned int mu,unsigned int i)
//{
//    unsigned int r;
//    r = 0;
//    for (int k = dim-1; k>=0; k--){
//        if( bit(mu, k) ) r = (r << 1) | bit(i,k);
//    }
//    return r;
//}



void Patch::initExchParticles(SmileiMPI* smpi, int ispec, PicParams& params, int tnum, int iDimOld)
{
    SmileiMPI_Cart2D* smpi2D = static_cast<SmileiMPI_Cart2D*>(smpi);
    Particles &cuParticles = (*vecSpecies[ispec]->particles);
    
    std::vector< std::vector<int> >* indexes_of_particles_to_exchange_per_thd = &vecSpecies[ispec]->indexes_of_particles_to_exchange_per_thd;
    std::vector<int>*                indexes_of_particles_to_exchange         = &vecSpecies[ispec]->indexes_of_particles_to_exchange;
    
    /********************************************************************************/
    // Build lists of indexes of particle to exchange per neighbor
    // Computed from indexes_of_particles_to_exchange computed during particles' BC
    /********************************************************************************/
    (*indexes_of_particles_to_exchange).clear();
        
    int tmp = 0;
    for (int tid=0 ; tid < (int)indexes_of_particles_to_exchange_per_thd->size() ; tid++)
	tmp += ((*indexes_of_particles_to_exchange_per_thd)[tid]).size();
    (*indexes_of_particles_to_exchange).resize( tmp );
        
    int k=0;
    for (int tid=0 ; tid < (int)indexes_of_particles_to_exchange_per_thd->size() ; tid++) {
	for (int ipart = 0 ; ipart < (int) ((*indexes_of_particles_to_exchange_per_thd)[tid]).size() ; ipart++ ) {
	    (*indexes_of_particles_to_exchange)[k] =  (*indexes_of_particles_to_exchange_per_thd)[tid][ipart] ;
	    k++;
	}
	((*indexes_of_particles_to_exchange_per_thd))[tid].clear();
    }
    sort( (*indexes_of_particles_to_exchange).begin(), (*indexes_of_particles_to_exchange).end() );
        
    int n_part_send = (*indexes_of_particles_to_exchange).size();
    int n_part_recv;
        
    int iPart;
    int n_particles;

    // Define where particles are going 
    for (int i=0 ; i<n_part_send ; i++) {
	iPart = (*indexes_of_particles_to_exchange)[i];

	if ( cuParticles.position(0,iPart) < min_local[0]) { 
	    if ( cuParticles.position(1,iPart) < min_local[1]) {
		vecSpecies[ispec]->specMPI.corner_buff_index_send[0][0].push_back( iPart );
	    }
	    else if ( cuParticles.position(1,iPart) >= max_local[1]) {
		vecSpecies[ispec]->specMPI.corner_buff_index_send[0][1].push_back( iPart );
	    }
	    else {
		vecSpecies[ispec]->specMPI.patch_buff_index_send[0][0].push_back( iPart );
	    }
	}
	else if ( cuParticles.position(0,iPart) >= max_local[0]) { 
	    if ( cuParticles.position(1,iPart) < min_local[1]) {
		vecSpecies[ispec]->specMPI.corner_buff_index_send[1][0].push_back( iPart );
	    }
	    else if ( cuParticles.position(1,iPart) >= max_local[1]) {
		vecSpecies[ispec]->specMPI.corner_buff_index_send[1][1].push_back( iPart );
	    }
	    else {
		vecSpecies[ispec]->specMPI.patch_buff_index_send[0][1].push_back( iPart );
	    }
	}
	else {
	    if ( cuParticles.position(1,iPart) < min_local[1]) {
		vecSpecies[ispec]->specMPI.patch_buff_index_send[1][0].push_back( iPart );
	    }
	    else if ( cuParticles.position(1,iPart) >= max_local[1]) {
		vecSpecies[ispec]->specMPI.patch_buff_index_send[1][1].push_back( iPart );
	    }
	    else {
		//If partciles is in but here, to be suppressed (= supp BC)
	    }
	}
    }
        
    /********************************************************************************/
    // Exchange number of particles to exchange to establish or not a communication
    /********************************************************************************/
    for (int iDim=0 ; iDim<2 ; iDim++) {
	for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
	    if (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {
		//n_part_send = (vecSpecies[ispec]->specMPI.patch_buff_index_send[iDim][iNeighbor]).size();
		vecSpecies[ispec]->specMPI.patch_buff_index_send_sz[iDim][iNeighbor] = (vecSpecies[ispec]->specMPI.patch_buff_index_send[iDim][iNeighbor]).size();
		int tag = buildtag( hindex, neighbor_[iDim][iNeighbor]);
		//MPI_Isend( &(vecSpecies[ispec]->specMPI.patch_buff_index_send_sz[iDim][iNeighbor]), 1, MPI_INT, 0, tag, MPI_COMM_SELF, &(vecSpecies[ispec]->specMPI.patch_srequest[iDim][iNeighbor]) );
		MPI_Isend( &(vecSpecies[ispec]->specMPI.patch_buff_index_send_sz[iDim][iNeighbor]), 1, MPI_INT, MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &(vecSpecies[ispec]->specMPI.patch_srequest[iDim][iNeighbor]) );
		//cout << hindex << " will sent " << n_part_send << " to " << neighbor_[iDim][iNeighbor] << " with tag " << tag << endl;
	    } // END of Send
	    else
		n_part_send = 0;
	    if (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
		vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[iDim][(iNeighbor+1)%2] = 0;
		int tag = buildtag( neighbor_[iDim][(iNeighbor+1)%2], hindex);
		//MPI_Irecv( &(vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[iDim][(iNeighbor+1)%2]), 1, MPI_INT, 0, tag, MPI_COMM_SELF, &(vecSpecies[ispec]->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]) );
		MPI_Irecv( &(vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[iDim][(iNeighbor+1)%2]), 1, MPI_INT, MPI_neighbor_[iDim][(iNeighbor+1)%2], tag, MPI_COMM_WORLD, &(vecSpecies[ispec]->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]) );
		//cout << hindex << " will recv from " << neighbor_[iDim][(iNeighbor+1)%2] << " with tag " << tag << endl;
	    }
	    else 
		vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[iDim][(iNeighbor+1)%2] = 0;
	}
    }

    for (int iDim=0 ; iDim<2 ; iDim++) {
	for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
	    if (corner_neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {
		//n_part_send = (vecSpecies[ispec]->specMPI.corner_buff_index_send[iDim][iNeighbor]).size();
		vecSpecies[ispec]->specMPI.corner_buff_index_send_sz[iDim][iNeighbor] = (vecSpecies[ispec]->specMPI.corner_buff_index_send[iDim][iNeighbor]).size();
		int tag = buildtag( hindex, corner_neighbor_[iDim][iNeighbor]);
		//MPI_Isend( &(vecSpecies[ispec]->specMPI.corner_buff_index_send_sz[iDim][iNeighbor]), 1, MPI_INT, 0, tag, MPI_COMM_SELF, &(vecSpecies[ispec]->specMPI.corner_srequest[iDim][iNeighbor]) );
		MPI_Isend( &(vecSpecies[ispec]->specMPI.corner_buff_index_send_sz[iDim][iNeighbor]), 1, MPI_INT, MPI_corner_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &(vecSpecies[ispec]->specMPI.corner_srequest[iDim][iNeighbor]) );
		//cout << hindex << " will sent " << n_part_send << " to " << corner_neighbor_[iDim][iNeighbor] << " with tag " << tag << endl;
	    } // END of Send
	    else
		n_part_send = 0;
	    if (corner_neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
		vecSpecies[ispec]->specMPI.corner_buff_index_recv_sz[iDim][(iNeighbor+1)%2] = 0;
		int tag = buildtag( corner_neighbor_[iDim][(iNeighbor+1)%2], hindex);
		//MPI_Irecv( &(vecSpecies[ispec]->specMPI.corner_buff_index_recv_sz[iDim][(iNeighbor+1)%2]), 1, MPI_INT, 0, tag, MPI_COMM_SELF, &(vecSpecies[ispec]->specMPI.corner_rrequest[iDim][(iNeighbor+1)%2]) );
		MPI_Irecv( &(vecSpecies[ispec]->specMPI.corner_buff_index_recv_sz[iDim][(iNeighbor+1)%2]), 1, MPI_INT, MPI_corner_neighbor_[iDim][(iNeighbor+1)%2], tag, MPI_COMM_WORLD, &(vecSpecies[ispec]->specMPI.corner_rrequest[iDim][(iNeighbor+1)%2]) );
		//cout << hindex << " will recv from " << corner_neighbor_[iDim][(iNeighbor+1)%2] << " with tag " << tag << endl;
	    }
	    else 
		vecSpecies[ispec]->specMPI.corner_buff_index_recv_sz[iDim][(iNeighbor+1)%2] = 0;
	}
    }

}


void Patch::initCommParticles(SmileiMPI* smpi, int ispec, PicParams& params, int tnum, int iDimOld)
{
    SmileiMPI_Cart2D* smpi2D = static_cast<SmileiMPI_Cart2D*>(smpi);
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
		    if (smpi2D->periods_[iDim]==1) {
			// Enabled periodicity
			if ( ( iNeighbor==0 ) &&  (Pcoordinates[iDim] == 0 ) &&( cuParticles.position(iDim,vecSpecies[ispec]->specMPI.patch_buff_index_send[iDim][iNeighbor][iPart]) < 0. ) ) {
			    cuParticles.position(iDim,vecSpecies[ispec]->specMPI.patch_buff_index_send[iDim][iNeighbor][iPart])     += x_max;
			}
			else if ( ( iNeighbor==1 ) &&  (Pcoordinates[iDim] == smpi2D->number_of_procs[iDim]-1 ) && ( cuParticles.position(iDim,vecSpecies[ispec]->specMPI.patch_buff_index_send[iDim][iNeighbor][iPart]) >= x_max ) ) {
			    cuParticles.position(iDim,vecSpecies[ispec]->specMPI.patch_buff_index_send[iDim][iNeighbor][iPart])     -= x_max;
			}
		    }
		    cuParticles.cp_particle(vecSpecies[ispec]->specMPI.patch_buff_index_send[iDim][iNeighbor][iPart], vecSpecies[ispec]->specMPI.patchVectorSend[iDim][iNeighbor]);
		}

		int tag = buildtag( hindex, neighbor_[iDim][iNeighbor]);
		typePartSend = smpi2D->createMPIparticles( &(vecSpecies[ispec]->specMPI.patchVectorSend[iDim][iNeighbor]), nbrOfProp );
		//MPI_Isend( &((vecSpecies[ispec]->specMPI.patchVectorSend[iDim][iNeighbor]).position(0,0)), 1, typePartSend, 0, tag, MPI_COMM_SELF, &(vecSpecies[ispec]->specMPI.patch_srequest[iDim][iNeighbor]) );
		MPI_Isend( &((vecSpecies[ispec]->specMPI.patchVectorSend[iDim][iNeighbor]).position(0,0)), 1, typePartSend, MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &(vecSpecies[ispec]->specMPI.patch_srequest[iDim][iNeighbor]) );
		//cout << hindex << " really send " << n_part_send << " to " << neighbor_[iDim][iNeighbor] << " with tag " << tag << endl;
		MPI_Type_free( &typePartSend );

	    } // END of Send
                
	    n_part_recv = vecSpecies[ispec]->specMPI.patch_buff_index_recv_sz[iDim][(iNeighbor+1)%2];
	    if ( (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
		typePartRecv = smpi2D->createMPIparticles( &(vecSpecies[ispec]->specMPI.patchVectorRecv[iDim][(iNeighbor+1)%2]), nbrOfProp );
		int tag = buildtag( neighbor_[iDim][(iNeighbor+1)%2], hindex);
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
		    if (smpi2D->periods_[iDim]==1) {
			// Enabled periodicity
			if ( (Pcoordinates[(iDim+1)%2] == 0 ) &&( cuParticles.position((iDim+1)%2,vecSpecies[ispec]->specMPI.corner_buff_index_send[iDim][iNeighbor][iPart]) < 0. ) ) {
			    cuParticles.position((iDim+1)%2,vecSpecies[ispec]->specMPI.corner_buff_index_send[iDim][iNeighbor][iPart])     += x_max;
			}
			if ( (Pcoordinates[iDim] == 0 ) &&( cuParticles.position(iDim,vecSpecies[ispec]->specMPI.corner_buff_index_send[iDim][iNeighbor][iPart]) < 0. ) ) {
			    cuParticles.position(iDim,vecSpecies[ispec]->specMPI.corner_buff_index_send[iDim][iNeighbor][iPart])     += y_max;
			}
			if ( (Pcoordinates[(iDim+1)%2] == smpi2D->number_of_procs[(iDim+1)%2]-1 ) && ( cuParticles.position((iDim+1)%2,vecSpecies[ispec]->specMPI.corner_buff_index_send[iDim][iNeighbor][iPart]) >= x_max ) ) {
			    cuParticles.position((iDim+1)%2,vecSpecies[ispec]->specMPI.corner_buff_index_send[iDim][iNeighbor][iPart])     -= x_max;
			}
			if ( (Pcoordinates[iDim] == smpi2D->number_of_procs[iDim]-1 ) && ( cuParticles.position(iDim,vecSpecies[ispec]->specMPI.corner_buff_index_send[iDim][iNeighbor][iPart]) >= y_max ) ) {
			    cuParticles.position(iDim,vecSpecies[ispec]->specMPI.corner_buff_index_send[iDim][iNeighbor][iPart])     -= y_max;
			}

		    }
		    cuParticles.cp_particle(vecSpecies[ispec]->specMPI.corner_buff_index_send[iDim][iNeighbor][iPart], vecSpecies[ispec]->specMPI.cornerVectorSend[iDim][iNeighbor]);
		}

		typePartSend = smpi2D->createMPIparticles( &(vecSpecies[ispec]->specMPI.cornerVectorSend[iDim][iNeighbor]), nbrOfProp );
		int tag = buildtag( hindex, corner_neighbor_[iDim][iNeighbor]);
		//MPI_Isend( &((vecSpecies[ispec]->specMPI.cornerVectorSend[iDim][iNeighbor]).position(0,0)), 1, typePartSend, 0, tag, MPI_COMM_SELF, &(vecSpecies[ispec]->specMPI.corner_srequest[iDim][iNeighbor]) );
		MPI_Isend( &((vecSpecies[ispec]->specMPI.cornerVectorSend[iDim][iNeighbor]).position(0,0)), 1, typePartSend, MPI_corner_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &(vecSpecies[ispec]->specMPI.corner_srequest[iDim][iNeighbor]) );
		//cout << hindex << " really send " << n_part_send << " to " << corner_neighbor_[iDim][iNeighbor] << " with tag " << tag << endl;
		MPI_Type_free( &typePartSend );

	    } // END of Send
                
	    n_part_recv = vecSpecies[ispec]->specMPI.corner_buff_index_recv_sz[iDim][(iNeighbor+1)%2];
	    if ( (corner_neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) && (n_part_recv!=0) ) {
		typePartRecv = smpi2D->createMPIparticles( &(vecSpecies[ispec]->specMPI.cornerVectorRecv[iDim][(iNeighbor+1)%2]), nbrOfProp );
		int tag = buildtag( corner_neighbor_[iDim][(iNeighbor+1)%2], hindex);
		//MPI_Irecv( &((vecSpecies[ispec]->specMPI.cornerVectorRecv[iDim][(iNeighbor+1)%2]).position(0,0)), 1, typePartRecv, 0, tag, MPI_COMM_SELF, &(vecSpecies[ispec]->specMPI.corner_rrequest[iDim][(iNeighbor+1)%2]) );
		MPI_Irecv( &((vecSpecies[ispec]->specMPI.cornerVectorRecv[iDim][(iNeighbor+1)%2]).position(0,0)), 1, typePartRecv, MPI_corner_neighbor_[iDim][(iNeighbor+1)%2], tag, MPI_COMM_WORLD, &(vecSpecies[ispec]->specMPI.corner_rrequest[iDim][(iNeighbor+1)%2]) );
		//cout << hindex << " will really recv from " << corner_neighbor_[iDim][(iNeighbor+1)%2] << " with tag " << tag << endl;
		MPI_Type_free( &typePartRecv );

	    } // END of Recv
                
	} // END for iNeighbor
    } // END for iDim


}

void Patch::finalizeCommParticles(SmileiMPI* smpi, int ispec, PicParams& params, int tnum, int iDimOld)
{
    SmileiMPI_Cart2D* smpi2D = static_cast<SmileiMPI_Cart2D*>(smpi);
    Particles &cuParticles = (*vecSpecies[ispec]->particles);


    std::vector< std::vector<int> >* indexes_of_particles_to_exchange_per_thd = &vecSpecies[ispec]->indexes_of_particles_to_exchange_per_thd;
    std::vector<int>*                indexes_of_particles_to_exchange         = &vecSpecies[ispec]->indexes_of_particles_to_exchange;

    std::vector<int>* cubmin = &vecSpecies[ispec]->bmin;
    std::vector<int>* cubmax = &vecSpecies[ispec]->bmax;

    int nmove,lmove; // local, OK
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
    int ii, iPart;
    // Push lost particles at the end of bins
    //! \todo For loop on bins, can use openMP here.
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
}


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


void Patch::initSumRhoJ( ElectroMagn* EMfields )
{
    // sum total charge density and currents
    initSumField( EMfields->rho_ );
    initSumField( EMfields->Jx_ );
    initSumField( EMfields->Jy_ );
    initSumField( EMfields->Jz_ );

}

void Patch::initSumField( Field* field )
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
    
    for (int iDim=0 ; iDim<patch_ndims_ ; iDim++) {
        for (int iNeighbor=0 ; iNeighbor<patch_nbNeighbors_ ; iNeighbor++) {
            std::vector<unsigned int> tmp(patch_ndims_,0);
            tmp[0] =    iDim  * n_elem[0] + (1-iDim) * oversize2[0];
            tmp[1] = (1-iDim) * n_elem[1] +    iDim  * oversize2[1];
            buf[iDim][iNeighbor].allocateDims( tmp );
        }
    }
     for (int iDim=0 ; iDim<patch_ndims_ ; iDim++) {
        for (int iNeighbor=0 ; iNeighbor<patch_nbNeighbors_ ; iNeighbor++) {
            std::vector<unsigned int> tmp(patch_ndims_,0);
            tmp[0] = 1 + 2 * patch_oversize[0] + isDual[0];
            tmp[1] = 1 + 2 * patch_oversize[1] + isDual[1];
            corner_buf[iDim][iNeighbor].allocateDims( tmp );
        }
    }
     
    int istart, ix, iy;
    /********************************************************************************/
    // Send/Recv in a buffer data to sum
    /********************************************************************************/
    for (int iDim=0 ; iDim<2 ; iDim++) {
        
	MPI_Datatype ntype = ntypeSum_[iDim][isDual[0]][isDual[1]];
	MPI_Request srequest[patch_ndims_][2];
	MPI_Request rrequest[patch_ndims_][2];
        
	for (int iNeighbor=0 ; iNeighbor<patch_nbNeighbors_ ; iNeighbor++) {
            
	    if (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {
		istart = iNeighbor * ( n_elem[iDim]- oversize2[iDim] ) + (1-iNeighbor) * ( 0 );
		ix = (1-iDim)*istart;
		iy =    iDim *istart;
		int tag = buildtag( hindex, neighbor_[iDim][iNeighbor]);
		//cout << hindex << " send to " << neighbor_[iDim][iNeighbor] << endl;
		//MPI_Isend( &(f2D->data_2D[ix][iy]), 1, ntype, 0, tag, MPI_COMM_SELF, &(f2D->specMPI.patch_srequest[iDim][iNeighbor]) );
		MPI_Isend( &(f2D->data_2D[ix][iy]), 1, ntype, MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &(f2D->specMPI.patch_srequest[iDim][iNeighbor]) );
	    } // END of Send
            
	    if (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
		int tmp_elem = (buf[iDim][(iNeighbor+1)%2]).dims_[0]*(buf[iDim][(iNeighbor+1)%2]).dims_[1];
		int tag = buildtag( neighbor_[iDim][(iNeighbor+1)%2], hindex);
		//cout << hindex << " recv from " << neighbor_[iDim][(iNeighbor+1)%2] << " ; n_elements = " << tmp_elem << endl;
		//MPI_Irecv( &( (buf[iDim][(iNeighbor+1)%2]).data_2D[0][0] ), tmp_elem, MPI_DOUBLE, 0, tag, MPI_COMM_SELF, &(f2D->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]) );
		MPI_Irecv( &( (buf[iDim][(iNeighbor+1)%2]).data_2D[0][0] ), tmp_elem, MPI_DOUBLE, MPI_neighbor_[iDim][(iNeighbor+1)%2], tag, MPI_COMM_WORLD, &(f2D->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]) );
	    } // END of Recv
            
	} // END for iNeighbor
    }

    for (int iDim=0 ; iDim<2 ; iDim++) {
	
	MPI_Datatype ntype = corner_ntypeSum_[0][isDual[0]][isDual[1]]; // 1st dimension useless
      
	for (int iNeighbor=0 ; iNeighbor<patch_nbNeighbors_ ; iNeighbor++) {
            
	    if (corner_neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {
		ix = iDim      * ( n_elem[0]     - oversize2[0]      ) + (1-iDim     ) * ( 0 );
		iy = iNeighbor * ( n_elem[1]- oversize2[1] ) + (1-iNeighbor) * ( 0 );
		int tag = buildtag( hindex, corner_neighbor_[iDim][iNeighbor]);
		int tabsize(0);
		MPI_Type_size( ntype, &tabsize );
		//cout << hindex << " send in diagonal to " << corner_neighbor_[iDim][iNeighbor] << " from " << ix << " " << iy << " " << tabsize <<  endl;
		//MPI_Isend( &(f2D->data_2D[ix][iy]), 1, ntype, 0, tag, MPI_COMM_SELF, &(f2D->specMPI.corner_srequest[iDim][iNeighbor]) );
		MPI_Isend( &(f2D->data_2D[ix][iy]), 1, ntype, MPI_corner_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &(f2D->specMPI.corner_srequest[iDim][iNeighbor]) );
	    } // END of Send
            
	    if (corner_neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
		int tmp_elem = (corner_buf[iDim][(iNeighbor+1)%2]).dims_[0]*(corner_buf[iDim][(iNeighbor+1)%2]).dims_[1];
		int tag = buildtag( corner_neighbor_[iDim][(iNeighbor+1)%2], hindex);
		//cout << hindex << " recv from " << corner_neighbor_[iDim][(iNeighbor+1)%2] << " ; n_elements = " << tmp_elem << endl;
		//MPI_Irecv( &( (corner_buf[iDim][(iNeighbor+1)%2]).data_2D[0][0] ), tmp_elem, MPI_DOUBLE, 0, tag, MPI_COMM_SELF, &(f2D->specMPI.corner_rrequest[iDim][(iNeighbor+1)%2]) );
		MPI_Irecv( &( (corner_buf[iDim][(iNeighbor+1)%2]).data_2D[0][0] ), tmp_elem, MPI_DOUBLE, MPI_corner_neighbor_[iDim][(iNeighbor+1)%2], tag, MPI_COMM_WORLD, &(f2D->specMPI.corner_rrequest[iDim][(iNeighbor+1)%2]) );
	    } // END of Recv
            
	} // END for iNeighbor
    }

} // END initSumField


void Patch::finalizeSumRhoJ( ElectroMagn* EMfields )
{
    // sum total charge density and currents
    finalizeSumField( EMfields->rho_ );
    finalizeSumField( EMfields->Jx_ );
    finalizeSumField( EMfields->Jy_ );
    finalizeSumField( EMfields->Jz_ );

}

void Patch::finalizeSumField( Field* field )
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

    for (int iDim=0 ; iDim<2 ; iDim++) {
 
        MPI_Datatype ntype = ntypeSum_[iDim][isDual[0]][isDual[1]];
        MPI_Status sstat    [patch_ndims_][2];
        MPI_Status rstat    [patch_ndims_][2];
	
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
            if (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {
		//cout << hindex << " is waiting for send at " << neighbor_[iDim][iNeighbor] << endl;
                MPI_Wait( &(f2D->specMPI.patch_srequest[iDim][iNeighbor]), &(sstat[iDim][iNeighbor]) );
            }
            if (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
		//cout << hindex << " is waiting for recv from " << neighbor_[iDim][(iNeighbor+1)%2] << endl;	
                MPI_Wait( &(f2D->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]), &(rstat[iDim][(iNeighbor+1)%2]) );
            }
        }
    }

    for (int iDim=0 ; iDim<2 ; iDim++) {
 
        MPI_Datatype ntype = corner_ntypeSum_[0][isDual[0]][isDual[1]];; // 1st dimension useless
        MPI_Status corner_sstat    [patch_ndims_][2];
        MPI_Status corner_rstat    [patch_ndims_][2];
	
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
            if (corner_neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) {
		//cout << hindex << " is waiting for corner send at " << corner_neighbor_[iDim][iNeighbor] << endl;
                MPI_Wait( &(f2D->specMPI.corner_srequest[iDim][iNeighbor]), &(corner_sstat[iDim][iNeighbor]) );
		//cout << hindex << " is waiting for corner send at " << corner_neighbor_[iDim][iNeighbor] << " ACHIEVED" << endl;
            }
            if (corner_neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
		//cout << hindex << " is waiting for corner recv from " << corner_neighbor_[iDim][(iNeighbor+1)%2] << endl;	
                MPI_Wait( &(f2D->specMPI.corner_rrequest[iDim][(iNeighbor+1)%2]), &(corner_rstat[iDim][(iNeighbor+1)%2]) );
            }
        }
    }

    /********************************************************************************/
    // Sum data on each process, same operation on both side
    /********************************************************************************/
    for (int iDim=0 ; iDim<2 ; iDim++) {
       
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
    }
        
    for (int iDim=0 ; iDim<2 ; iDim++) {
       
        for (int iNeighbor=0 ; iNeighbor<nbNeighbors_ ; iNeighbor++) {
	    int ix0 = iDim      * ( n_elem[0]     - oversize2[0]      ) + (1-iDim     ) * ( 0 );
	    int iy0 = iNeighbor * ( n_elem[1]- oversize2[1] ) + (1-iNeighbor) * ( 0 );
            if (corner_neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) {
                for (unsigned int ix=0 ; ix< (corner_buf[iDim][(iNeighbor+1)%2]).dims_[0] ; ix++) {
                    for (unsigned int iy=0 ; iy< (corner_buf[iDim][(iNeighbor+1)%2]).dims_[1] ; iy++)
                        f2D->data_2D[ix0+ix][iy0+iy] += (corner_buf[iDim][(iNeighbor+1)%2])(ix,iy);
                }
            } // END if
            
        } // END for iNeighbor
    }       

    for (int iDim=0 ; iDim<patch_ndims_ ; iDim++) {
        for (int iNeighbor=0 ; iNeighbor<patch_nbNeighbors_ ; iNeighbor++) {
	    buf[iDim][iNeighbor].deallocateDims();
        }
    }
     for (int iDim=0 ; iDim<patch_ndims_ ; iDim++) {
        for (int iNeighbor=0 ; iNeighbor<patch_nbNeighbors_ ; iNeighbor++) {
	    corner_buf[iDim][iNeighbor].deallocateDims();
        }
    }

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

            if ( (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) ) {

                istart = iNeighbor * ( n_elem[iDim]- (2*patch_oversize[iDim]+1+isDual[iDim]) ) + (1-iNeighbor) * ( 2*patch_oversize[iDim] + isDual[iDim] );
                ix = (1-iDim)*istart;
                iy =    iDim *istart;
		int tag = buildtag( hindex, neighbor_[iDim][iNeighbor]);
                //MPI_Isend( &(f2D->data_2D[ix][iy]), 1, ntype, 0, tag, MPI_COMM_SELF, &(f2D->specMPI.patch_srequest[iDim][iNeighbor]) );
                MPI_Isend( &(f2D->data_2D[ix][iy]), 1, ntype, MPI_neighbor_[iDim][iNeighbor], tag, MPI_COMM_WORLD, &(f2D->specMPI.patch_srequest[iDim][iNeighbor]) );

            } // END of Send

            if ( (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL) ) {

                istart = ( (iNeighbor+1)%2 ) * ( n_elem[iDim] - 1 ) + (1-(iNeighbor+1)%2) * ( 0 )  ;
                ix = (1-iDim)*istart;
                iy =    iDim *istart;
 		int tag = buildtag( neighbor_[iDim][(iNeighbor+1)%2], hindex);
		//MPI_Irecv( &(f2D->data_2D[ix][iy]), 1, ntype, 0, tag, MPI_COMM_SELF, &(f2D->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]));
		MPI_Irecv( &(f2D->data_2D[ix][iy]), 1, ntype, MPI_neighbor_[iDim][(iNeighbor+1)%2], tag, MPI_COMM_WORLD, &(f2D->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]));

            } // END of Recv

        } // END for iNeighbor

    } // END for iDim

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
            if ( (neighbor_[iDim][iNeighbor]!=MPI_PROC_NULL) ) {
                MPI_Wait( &(f2D->specMPI.patch_srequest[iDim][iNeighbor]), &(sstat[iDim][iNeighbor]) );
            }
            if ( (neighbor_[iDim][(iNeighbor+1)%2]!=MPI_PROC_NULL)  ) {
                MPI_Wait( &(f2D->specMPI.patch_rrequest[iDim][(iNeighbor+1)%2]), &(rstat[iDim][(iNeighbor+1)%2]) );
            }
        }

    } // END for iDim
}



VectorPatch::VectorPatch()
{
}

VectorPatch::~VectorPatch()
{
}

void VectorPatch::exchangeParticles(int ispec, PicParams &params, SmileiMPI* smpi)
{
    int useless(0);

    for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++) {
	(*this)(ipatch)->initExchParticles(smpi, ispec, params, useless, useless);
    }
    for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++) {
	(*this)(ipatch)->initCommParticles(smpi, ispec, params, useless, useless);
    }
    for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++) {
	(*this)(ipatch)->finalizeCommParticles(smpi, ispec, params, useless, useless);
	(*this)(ipatch)->vecSpecies[ispec]->sort_part();
    }

}

void VectorPatch::sumRhoJ( int ispec )
{

    for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++) {
	(*this)(ipatch)->initSumField( (*this)(ipatch)->EMfields->rho_ ); // initialize
    }
    for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++) {
	(*this)(ipatch)->finalizeSumField( (*this)(ipatch)->EMfields->rho_ ); // finalize (waitall + sum)
    }

    for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++) {
	(*this)(ipatch)->initSumField( (*this)(ipatch)->EMfields->Jx_ ); // initialize
    }
    for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++) {
	(*this)(ipatch)->finalizeSumField( (*this)(ipatch)->EMfields->Jx_ ); // finalize (waitall + sum)
    }
    for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++) {
	(*this)(ipatch)->initSumField( (*this)(ipatch)->EMfields->Jy_ ); // initialize
    }
    for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++) {
	(*this)(ipatch)->finalizeSumField( (*this)(ipatch)->EMfields->Jy_ ); // finalize (waitall + sum)
    }
    for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++) {
	(*this)(ipatch)->initSumField( (*this)(ipatch)->EMfields->Jz_ ); // initialize
    }
    for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++) {
	(*this)(ipatch)->finalizeSumField( (*this)(ipatch)->EMfields->Jz_ ); // finalize (waitall + sum)
    }


}

void VectorPatch::exchangeE( )
{
    for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++)
	(*this)(ipatch)->initExchange( (*this)(ipatch)->EMfields->Ex_ );
    for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++)
	(*this)(ipatch)->finalizeExchange( (*this)(ipatch)->EMfields->Ex_ );

    for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++)
	(*this)(ipatch)->initExchange( (*this)(ipatch)->EMfields->Ey_ );
    for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++)
	(*this)(ipatch)->finalizeExchange( (*this)(ipatch)->EMfields->Ey_ );

    for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++)
	(*this)(ipatch)->initExchange( (*this)(ipatch)->EMfields->Ez_ );
    for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++)
	(*this)(ipatch)->finalizeExchange( (*this)(ipatch)->EMfields->Ez_ );

}

void VectorPatch::exchangeB( )
{
    for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++)
	(*this)(ipatch)->initExchange( (*this)(ipatch)->EMfields->Bx_ );
    for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++)
	(*this)(ipatch)->finalizeExchange( (*this)(ipatch)->EMfields->Bx_ );

    for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++)
	(*this)(ipatch)->initExchange( (*this)(ipatch)->EMfields->By_ );
    for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++)
	(*this)(ipatch)->finalizeExchange( (*this)(ipatch)->EMfields->By_ );

    for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++)
	(*this)(ipatch)->initExchange( (*this)(ipatch)->EMfields->Bz_ );
    for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++)
	(*this)(ipatch)->finalizeExchange( (*this)(ipatch)->EMfields->Bz_ );

}

void VectorPatch::computeGlobalDiags(int timestep)
{
    
    computeScalarsDiags(timestep);
    //computeGlobalDiags(probes); // HDF5 write done per patch in DiagProbes::*
    //computeGlobalDiags(phases);
}

void VectorPatch::computeScalarsDiags(int timestep)
{
    int scalars_every( (*this)(0)->Diags->scalars.every );
    if (timestep % scalars_every != 0) return;

    int nDiags( (*this)(0)->Diags->scalars.out_list.size() );
    // Initialize scalars iterator on 1st diag
    for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++)
	(*this)(ipatch)->Diags->scalars.itDiagScalar =  (*this)(ipatch)->Diags->scalars.out_list.begin();


    for (int idiags = 0 ; idiags<nDiags ; idiags++) {
	string diagName( (*this)(0)->Diags->scalars.itDiagScalar->first );

	if ( ( diagName.find("Min") == std::string::npos ) && ( diagName.find("Max") == std::string::npos ) ) {
	    double sum(0.);
	    for (unsigned int ipatch=0 ; ipatch<this->size() ; ipatch++) {
		sum += (*this)(ipatch)->Diags->scalars.itDiagScalar->second;
		if (ipatch)
		    (*this)(ipatch)->Diags->scalars.itDiagScalar++;
	    }
	    (*this)(0)->Diags->scalars.itDiagScalar->second = sum;
	    (*this)(0)->Diags->scalars.itDiagScalar++;
	}
	else if ( diagName.find("MinCell") != std::string::npos ) {
	    vector<pair<string,double> >::iterator iterVal    = (*this)(0)->Diags->scalars.itDiagScalar-1;
	    vector<pair<string,double> >::iterator iterValRef = (*this)(0)->Diags->scalars.itDiagScalar-1;
	    double min( iterValRef->second );

	    for (unsigned int ipatch=1 ; ipatch<this->size() ; ipatch++) {
		if ((*this)(ipatch)->Diags->scalars.itDiagScalar->second < min) {
		    min = (*this)(ipatch)->Diags->scalars.itDiagScalar->second;
		    iterVal = (*this)(ipatch)->Diags->scalars.itDiagScalar-1;
		}
		if (ipatch)
		    (*this)(ipatch)->Diags->scalars.itDiagScalar++;
	    }
	    (*this)(0)->Diags->scalars.itDiagScalar->second = min;
	    iterValRef->second = iterVal->second;

	    (*this)(0)->Diags->scalars.itDiagScalar++;	    
	}
	else if ( diagName.find("MaxCell") != std::string::npos ) {
	    vector<pair<string,double> >::iterator iterVal    = (*this)(0)->Diags->scalars.itDiagScalar-1;
	    vector<pair<string,double> >::iterator iterValRef = (*this)(0)->Diags->scalars.itDiagScalar-1;
	    double max( iterValRef->second );

	    for (unsigned int ipatch=1 ; ipatch<this->size() ; ipatch++) {
		if ((*this)(ipatch)->Diags->scalars.itDiagScalar->second > max) {
		    max = (*this)(ipatch)->Diags->scalars.itDiagScalar->second;
		    iterVal = (*this)(ipatch)->Diags->scalars.itDiagScalar-1;
		}
		if (ipatch)
		    (*this)(ipatch)->Diags->scalars.itDiagScalar++;
	    }
	    (*this)(0)->Diags->scalars.itDiagScalar->second = max;
	    iterValRef->second = iterVal->second;

	    (*this)(0)->Diags->scalars.itDiagScalar++;	    
	}

	// Go to next diag
    }

    // After MPI sync
    //(*this)(0)->Diags->scalars.write(timestep);

}
