#include "Patch.h"
#include <iostream>

using namespace std;

//Functions for manipulating Hilbert curves in an arbitrary number of dimensions and arbitrary size.
//Taken from Chris Hamilton, Technical Report CS-2006-07, Faculty of computer Science, Halifax.

//!Get kth bit of i.
unsigned int bit( unsigned int i, unsigned int k )
{
    return ( i>>k )&1 ;
}
//!Set kth bit of i to value.
void setbit( unsigned int *i, unsigned int k, unsigned int value )
{
    *i = ( *i & ~( 1<<k ) ) | ( value<<k );
    return;
}

//Bitwise rotation operators.
//Works only for shift <= dim.
unsigned int rotl( unsigned int value, unsigned int shift, unsigned int dim )
{
    //Evaluate left side of the rotated value then right side of it. Finally retain only the first dim bits.
    return ( ( value << shift ) | ( value >> ( dim - shift ) ) ) & ( ( 1<<dim ) - 1 );
}
unsigned int rotr( unsigned int value, unsigned int shift, unsigned int dim )
{
    //Evaluate right side of the rotated value then left side of it. Finally retain only the first dim bits.
    return ( ( value >> shift ) | ( value << ( dim - shift ) ) ) & ( ( 1<<dim ) - 1 );
}

//Generates the binary reflected Gray Code
unsigned int gc( unsigned int i )
{
    return i^( i>>1 );
}
//Given a non-negative integer g, calculates the non-negative integer i such that gc(i)=g.
unsigned int gcinv( unsigned int g )
{
    unsigned int i, j;
    i=g;
    j=1;
    while( ( ( ( unsigned int )1 )<<j ) <= g ) {
        i = i ^ ( g >> j );
        j++;
    }
    return i;
}

// Tsb = trailing set bit. It is the number of trailing set bits in the binary representation of i.
// tsb is also the inter sub-hypercube directio, g(i).
unsigned int tsb( unsigned int i )
{
    unsigned int k;
    k = 0;
    while( i & 1 ) {
        i = i>>1;
        k++;
    }
    return k;
}
// Direction computes the sequence of intra sub-hypercube direction, d(i) for 0 <= i < 2^dim.
unsigned int direction( unsigned int i, unsigned int dim )
{
    if( i == 0 ) {
        return 0;
    } else if( i & 1 ) {
        return tsb( i )%dim;
    } else {
        return tsb( i-1 )%dim;
    }
}
// Entry computes the sequence of entry points, e(i) for 0 <= i < 2^dim.
unsigned int entry( unsigned int i )
{
    if( i == 0 ) {
        return 0;
    } else {
        return gc( 2*( ( i-1 )/2 ) );
    }
}
//!Ted is the transformation such that the gc ordering of sub-hypercubes in the Hilbert curve defined by e and d will map tot he standard binary reflected gc.
void ted( unsigned int e, unsigned int d, unsigned int *b, unsigned int dim )
{
    *b = rotr( *b ^ e, d+1, dim );
    return;
}
void tedinv( unsigned int e, unsigned int d, unsigned int *b, unsigned int dim )
{
    *b = rotl( *b, d+1, dim ) ^ e;
    return;
}
//!Hilbert index2D calculates the Hilbert index h of a patch of coordinates x,y for a simulation box with 2^m patches per side (2^(2*m) patches in total).
unsigned int hilbertindex( unsigned int m, unsigned int x, unsigned int y, unsigned int *einit, unsigned int *dinit )
{
    unsigned int e, d, h, l, w;
    e = *einit;
    d = *dinit;
    h = 0;
    for( int i = m-1; i>=0; i-- ) {
        l = bit( y, i )*2 + bit( x, i ); //ith bit of y at the leftmost position of l, and ith bit of x at the rightmost position of l.
        ted( e, d, &l, 2 );
        w = gcinv( l );
        e = e ^ ( rotl( entry( w ), d+1, 2 ) );
        d = ( d + direction( w, 2 ) +1 )%2 ;
        h = ( h<<2 )|w;
    }
    
    *einit = e;
    *dinit = d;
    
    return h;
}
//!Hilbert index3D calculates the Hilbert index h of a patch of coordinates x,y,z for a simulation box with 2^m patches per side (2^(3*m) patches in total).
unsigned int hilbertindex( unsigned int m, unsigned int x, unsigned int y, unsigned int z, unsigned int einit, unsigned int dinit )
{
    unsigned int e, d, h, l, w;
    h = 0;
    e = einit;
    d = dinit;
    for( int i = m-1; i>=0; i-- ) {
        l = bit( z, i )*4 + bit( y, i )*2 + bit( x, i );
        ted( e, d, &l, 3 );
        w = gcinv( l );
        e = e ^ ( rotl( entry( w ), d+1, 3 ) );
        d = ( d + direction( w, 3 ) + 1 )%3 ;
        h = ( h<<3 )|w;
    }
    return h;
}

//!Hilbert index2D inv  calculates the coordinates x,y of the patch of Hilbert index h in a simulation box with 2^m patches per side (2^(2*m) patches in total).
void hilbertindexinv( unsigned int m, unsigned int *x, unsigned int *y, unsigned int h, unsigned int einit, unsigned int dinit )
{
    unsigned int e, d, l, w;
    e = einit;
    d = dinit;
    *x=0;
    *y=0;
    for( int i = m-1; i>=0 ; i-- ) {
        w = ( ( bit( h, 2*i+1 ) )<<1 ) + bit( h, 2*i );
        l = gc( w );
        tedinv( e, d, &l, 2 );
        setbit( x, ( unsigned int )i, bit( l, 0 ) );
        setbit( y, ( unsigned int )i, bit( l, 1 ) );
        e = e ^ ( rotl( entry( w ), d+1, 2 ) );
        d = ( d + direction( w, 2 ) +1 )%2 ;
    }
    return;
}
void hilbertindexinv( unsigned int m, unsigned int *x, unsigned int *y, unsigned int *z, unsigned int h, unsigned int einit, unsigned int dinit )
{
    unsigned int e, d, l, w;
    e = einit;
    d = dinit;
    *x=0;
    *y=0;
    *z=0;
    for( int i = m-1; i>=0 ; i-- ) {
        w = ( ( bit( h, 3*i+2 ) )<<2 ) + ( ( bit( h, 3*i+1 ) )<<1 ) + bit( h, 3*i );
        l = gc( w );
        tedinv( e, d, &l, 3 );
        setbit( x, ( unsigned int )i, bit( l, 0 ) );
        setbit( y, ( unsigned int )i, bit( l, 1 ) );
        setbit( z, ( unsigned int )i, bit( l, 2 ) );
        e = e ^ ( rotl( entry( w ), d+1, 3 ) );
        d = ( d + direction( w, 3 ) +1 )%3 ;
    }
    return;
}


//The "general" versions of the functions allow a different number of patch in each direction.

//!General Hilbert index2D calculates the  Hilbert index h of a patch of coordinates x,y for a simulation box with 2^mi patches per side (2^(m0+m1)) patches in total).
unsigned int generalhilbertindex( unsigned int m0, unsigned int m1, int x, int y, unsigned int *einit, unsigned int *dinit )
{

    if( ( x<0 ) || ( x>=( 1<<m0 ) ) || ( y<0 ) || ( y>=( 1<<m1 ) ) ) {
        return MPI_PROC_NULL ;
    }
    
    unsigned int h, mmin, mmax, l, localx, localy, *target;
    h=0;
    *dinit=0;
    localx = ( unsigned int )x;
    localy = ( unsigned int )y;
    if( m0 >= m1 ) {
        target = &localx;
        mmin = m1;
        mmax = m0;
    } else {
        target = &localy;
        *dinit=1;
        mmin = m0;
        mmax = m1;
    }
    for( int i= ( int )mmax-1; i >= ( int )mmin ; i-- ) {
        l = bit( *target, i );
        h += l*( 1<<( i+mmin ) );
        *target -= l*( 1<<i );
    }
    if( mmin > 0 ) {
        h += hilbertindex( mmin, localx, localy, einit, dinit );
    }
    return h;
}
unsigned int generalhilbertindex( unsigned int m0, unsigned int m1, int x, int y )
{

    if( ( x<0 ) || ( x>=( 1<<m0 ) ) || ( y<0 ) || ( y>=( 1<<m1 ) ) ) {
        return MPI_PROC_NULL ;
    }
    
    unsigned int h, l, localx, localy, *target, einit, dinit;
    int mmin, mmax;
    h=0;
    dinit=0;
    einit=0;
    localx = ( unsigned int )x;
    localy = ( unsigned int )y;
    if( m0 >= m1 ) {
        target = &localx;
        mmin = m1;
        mmax = m0;
    } else {
        target = &localy;
        dinit=1;
        mmin = m0;
        mmax = m1;
    }
    for( int i= mmax-1; i >= mmin ; i-- ) {
        l = bit( *target, i );
        h += l*( 1<<( i+mmin ) );
        *target -= l*( 1<<i );
    }
    if( mmin > 0 ) {
        h += hilbertindex( ( unsigned int )mmin, localx, localy, &einit, &dinit );
    }
    return h;
}

//!General Hilbert index3D calculates the compact Hilbert index h of a patch of coordinates x,y,z for a simulation box with 2^mi patches per side (2^(m0+m1+m2)) patches in total).
unsigned int generalhilbertindex( unsigned int m0, unsigned int m1, unsigned int m2,  int x,  int y,  int z )
{
    if( ( x<0 ) || ( x>=( 1<<m0 ) ) || ( y<0 ) || ( y>=( 1<<m1 ) ) || ( z<0 ) || ( z>=( 1<<m2 ) ) ) {
        return MPI_PROC_NULL ;
    }
    
    unsigned int h, e, d, *einit, *dinit, dimmin, dimmax, dimmed, mi[3], localp[3], tempp[3];
    h=0;
    e=0;
    d=0;
    dinit=&d;
    einit=&e;
    //Store positions and dimensions in arrays
    localp[0] = ( unsigned int )x;
    localp[1] = ( unsigned int )y;
    localp[2] = ( unsigned int )z;
    mi[0] = m0;
    mi[1] = m1;
    mi[2] = m2;
    //Compare dimension sizes
    if( ( m0 >= m1 ) && ( m0 >= m2 ) ) {
        dimmax = 0;
    } else if( ( m1 > m0 ) && ( m1 >= m2 ) ) {
        dimmax = 1;
    } else {
        dimmax = 2;
    }
    if( mi[( dimmax+1 )%3] >= mi[( dimmax+2 )%3] ) {
        dimmed = ( dimmax+1 )%3;
        dimmin = ( dimmax+2 )%3;
    } else {
        dimmed = ( dimmax+2 )%3;
        dimmin = ( dimmax+1 )%3;
    }
    // First approach on a flattened 2D grid along dimmax and dimmed. The 3D grid is projected along dimmin axis.
    tempp[dimmax] = localp[dimmax] >> mi[dimmin]; //Erase last mi[dimmin] bits. Not relevent for this phase.
    tempp[dimmed] = localp[dimmed] >> mi[dimmin]; //Erase last mi[dimmin] bits. Not relevent for this phase.
    h += generalhilbertindex( mi[dimmax]-mi[dimmin], mi[dimmed]-mi[dimmin], tempp[dimmax], tempp[dimmed], einit, dinit )*( 1<<( 3*mi[dimmin] ) );
    
    //Now in a local cube of side mi[dimmin]. The local entry point "einit" and initial direction "dinit" of the local hilbert curve has been
    //determined by the previous call to compacthilbertindex2.
    //Relative position in the local cube is given by the last mi[dimmin] bits of the position.
    tempp[dimmax] = localp[dimmax] & ( ( 1<<mi[dimmin] )-1 ); //Only keep the last mi[dimmin] bits.
    tempp[dimmed] = localp[dimmed] & ( ( 1<<mi[dimmin] )-1 ); //Only keep the last mi[dimmin] bits.
    tempp[dimmin] = localp[dimmin] & ( ( 1<<mi[dimmin] )-1 ); //Only keep the last mi[dimmin] bits.
    //Add local index to the previously calculated one.
    h += hilbertindex( mi[dimmin], tempp[dimmax], tempp[dimmed], tempp[dimmin], *einit, *dinit );
    
    
    
    return h;
}
//!General Hilbert index inverse calculates the coordinates x,y of a patch for a given Hilbert index h in a simulation box with 2^mi patches per side (2^(m0+m1) patches in total)
//2D version
void generalhilbertindexinv( unsigned int m0, unsigned int m1, unsigned int *x, unsigned int *y, unsigned int h )
{
    unsigned int einit, dinit, l, localh, *target, shift ;
    int mmin, mmax;
    einit = 0;
    dinit = 0;
    shift = 0;
    localh = h;
    //Compare dimensions. Target points at the dimension which must be shifted.
    if( m0 >= m1 ) {
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
    for( int i= mmax+mmin-1; i >= mmin+mmin ; i-- ) {
        l = bit( localh, i );
        shift += l*( 1<<( i-mmin ) );
        localh -= l*( 1<<i );
    }
    //Run the cubic inversion algorithm in the sub hypercube.
    hilbertindexinv( mmin, x, y, localh, einit, dinit );
    //Shift the appropriate coordinate by the necessary value.
    *target += shift;
    
    return;
    
}
//3D version
void generalhilbertindexinv( unsigned int m0, unsigned int m1, unsigned int m2,  unsigned int *x, unsigned int *y, unsigned int *z, unsigned int h )
{
    unsigned int e, d, dimmin, dimmax, dimmed, mi[3], *localp[3], tempp[3], localh;
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
    if( ( m0 >= m1 ) && ( m0 >= m2 ) ) {
        dimmax = 0;
    } else if( ( m1 > m0 ) && ( m1 >= m2 ) ) {
        dimmax = 1;
    } else {
        dimmax = 2;
    }
    if( mi[( dimmax+1 )%3] >= mi[( dimmax+2 )%3] ) {
        dimmed = ( dimmax+1 )%3;
        dimmin = ( dimmax+2 )%3;
    } else {
        dimmed = ( dimmax+2 )%3;
        dimmin = ( dimmax+1 )%3;
    }
    
    
    //Localize in which sub hypercube the point is. Do not account for the first 3*dimmin bits of h.
    localh = ( h >> ( mi[dimmin]*3 ) );
    //Run the 2D inversion algorithm on the reduced domain.
    generalhilbertindexinv( mi[dimmax]-mi[dimmin], mi[dimmed]-mi[dimmin], localp[dimmax], localp[dimmed], localh );
    // Now local P stores the position of the cube in the 2D domain
    // We need to run the 3D inversion algorithm on this cube with the correct entry point and direction.
    // Run the 2D indexgenerator in order to evaluate e and d.
    localh = generalhilbertindex( mi[dimmax]-mi[dimmin], mi[dimmed]-mi[dimmin], *localp[dimmax], *localp[dimmed], &e, &d );
    //Transform coordinates in the global frame.
    *localp[dimmax] *= ( 1<<mi[dimmin] );
    *localp[dimmed] *= ( 1<<mi[dimmin] );
    *localp[dimmin] = 0 ;
    //Use only first bits of h for the local hypercube.
    localh = h & ( ( 1<<( mi[dimmin]*3 ) )-1 );
    //Run the cubic inversion algorithm in the local sub hypercube.
    hilbertindexinv( mi[dimmin], &tempp[dimmax], &tempp[dimmed], &tempp[dimmin], localh, e, d );
    //Add results to the coordinates.
    *localp[dimmax] += tempp[dimmax];
    *localp[dimmed] += tempp[dimmed];
    *localp[dimmin] += tempp[dimmin];
    
    return;
}


