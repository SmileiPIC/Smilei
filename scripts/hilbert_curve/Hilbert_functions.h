#include <vector>


//!Get bit.
unsigned int bit(unsigned int i, unsigned int k);
//!Set bit.
void setbit(unsigned int* i, unsigned int k, unsigned int value);
//!Bitwise rotation operators.
unsigned int rotl(unsigned int value, unsigned int shift, unsigned int dim);
unsigned int rotr(unsigned int value, unsigned int shift, unsigned int dim);
//!Generates binary reflected Gray Code
unsigned int gc(unsigned int i);
//!Gray code inverse
unsigned int gcinv(unsigned int g);
//!Tsb = trailing set bit. It is the number of trailing set bits in the binary representation of i.
//!tsb is also the inter sub-hypercube directio, g(i).
unsigned int tsb(unsigned int i);
//!Direction computes the sequence of intra sub-hypercube direction, d(i) for 0 <= i < 2^dim.
unsigned int direction(unsigned int i, unsigned int dim);
//!Entry computes the sequence of entry points, e(i) for 0 <= i < 2^dim.
unsigned int entry(unsigned int i);
//!Ted is the transformation such that the gc ordering of sub-hypercubes in the Hilbert curve defined by e and d will map tot he standard binary reflected gc.
void ted(unsigned int e, unsigned int d, unsigned int *b, unsigned int dim);
void tedinv(unsigned int e, unsigned int d, unsigned int *b, unsigned int dim);
//!Hilbert index calculates the Hilbert index h of a patch of coordinates x,y(z) for a simulation box with 2^m patches per side (2^(2 or 3*m) patches in total).
unsigned int hilbertindex(unsigned int m, unsigned int x, unsigned int y, unsigned int *einit, unsigned int *dinit);
unsigned int hilbertindex(unsigned int m, unsigned int x, unsigned int y, unsigned int z, unsigned int einit, unsigned int dinit);
//!Hilbert index inv returns the coordinates of the patch of Hilbert index h in domain of size 2^(m*dim).
void hilbertindexinv(unsigned int m, unsigned int* x, unsigned int* y, unsigned int h, unsigned int einit, unsigned int dinit);
void hilbertindexinv(unsigned int m, unsigned int* x, unsigned int* y, unsigned int* z, unsigned int h, unsigned int einit, unsigned int dinit);
//!General Hilbert index returns the Hilbert index h of a patch of coordinates x,y,z for a simulation box with 2^mi patches per side (2^(m0+m1+m2) patches in total).
//The 2D version of this function stores the final entry point and direction in enit and dinit (needed by the 3D version).
unsigned int generalhilbertindex(unsigned int m0, unsigned int m1,  int x,  int y, unsigned int *einit, unsigned int *dinit);
unsigned int generalhilbertindex(unsigned int m0, unsigned int m1,  int x,  int y);
unsigned int generalhilbertindex(unsigned int m0, unsigned int m1, unsigned int m2,  int x,  int y,  int z);
//!General Hilbert index inv calculates the coordinates x,y,z of a patch for a given Hilbert index h in a simulation box with 2^mi patches per side (2^(m0+m1+m2) patches in total)
void generalhilbertindexinv(unsigned int m0, unsigned int m1, unsigned int* x, unsigned int* y, unsigned int h);
void generalhilbertindexinv(unsigned int m0, unsigned int m1, unsigned int m2, unsigned int* x, unsigned int* y, unsigned int* z, unsigned int h);


