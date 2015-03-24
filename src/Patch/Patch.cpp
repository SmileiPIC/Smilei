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

	Pcoordinates[0] = ipatch%m0;
	Pcoordinates[1] = ipatch/m0 ;
	std::cout << ipatch << " " << Pcoordinates[0] << " " << Pcoordinates[1] << std::endl;

	min_local.resize(params.nDim_field, 0.);
	max_local.resize(params.nDim_field, 0.);
	cell_starting_global_index.resize(params.nDim_field, 0);
	for (int i = 0 ; i<params.nDim_field ; i++) {
	    min_local[i] = smpi->getDomainLocalMin(i) + Pcoordinates[i]*params.n_space[i]*params.cell_length[i];
	    max_local[i] = min_local[i] + params.n_space[i]*params.cell_length[i];
	    cell_starting_global_index[i] += Pcoordinates[i]*params.n_space[i];
	}

	std::cout << "Create patch\n";


	vecSpecies = SpeciesFactory::createVector(params, smpi, this);

	/* // + min_loc/cell_index(ref smpi,  & sort) // OK through this 
 * 	   std::cout << "Patch created\n";*/
	// + new n_space -> in PatchFactory
	// patchID : ok through coord
	// create Pos : OK
	//return;

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
