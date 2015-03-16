#ifndef PATCH_H
#define PATCH_H

#include <vector>

#include "SpeciesFactory.h"
#include "ElectroMagnFactory.h"
#include "InterpolatorFactory.h"
#include "ProjectorFactory.h"

#include "PicParams.h"
#include "LaserParams.h"
#include "SmileiMPI.h"
#include <iostream>
#include <cstdlib>
#include <iomanip>


//! Class Patch : sub MPI domain 
//!     Collection of patch = MPI domain
class Patch
{

public:
    //! Constructor for Patch
    Patch(PicParams& params, LaserParams& laser_params, SmileiMPI* smpi, unsigned int m0, unsigned int m1, unsigned int m2, unsigned int ipatch) {

	vecSpecies = SpeciesFactory::createVector(params, smpi);              // + patchId + min_loc/cell_index(ref smpi, creta Pos & sort) + new n_space
	// -> partBoundCond : min/max_loc (smpi)
	EMfields   = ElectroMagnFactory::create(params, laser_params, smpi);  // + patchId + new n_space (now = params by smpi) + BC
	Interp     = InterpolatorFactory::create(params, smpi);               // + patchId -> idx_domain_begin (now = ref smpi)
	Proj       = ProjectorFactory::create(params, smpi);                  // + patchId -> idx_domain_begin (now = ref smpi)

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
	
    };

    //! Destructor for Patch
    ~Patch() {

	delete Proj;
	delete Interp;
	delete EMfields;
	for (unsigned int ispec=0 ; ispec<vecSpecies.size(); ispec++) delete vecSpecies[ispec];
	vecSpecies.clear();
	    
    };

    std::vector<Species*> vecSpecies;
    ElectroMagn* EMfields;

    Interpolator* Interp;
    Projector* Proj;

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
    unsigned int hilbertindex(unsigned int m, unsigned int x, unsigned int y);
    unsigned int hilbertindex(unsigned int m, unsigned int x, unsigned int y, unsigned int z);
    void hilbertindexinv(unsigned int m, unsigned int* x, unsigned int* y, unsigned int h);
    //!extractMask extracts a mask µ indicating which axes are active at a given iteration i of the compact hilbert index.
    unsigned int extractmask(unsigned int m0,unsigned int  m1, unsigned int i );
    unsigned int extractmask(unsigned int m0,unsigned int  m1, unsigned int  m2, unsigned int i );
    //!Gray Code Rank.
    unsigned int gcr(unsigned int dim, unsigned int mu,unsigned int i);
    //!Gray Code Rank Inverse.
    unsigned int gcrinv(unsigned int dim, unsigned int mu,unsigned int pi, unsigned int r);
    //!Hilbert index returns the Hilbert index h of a patch of coordinates x,y,z for a simulation box with 2^mi patches per side (2^(m0+m1+m2) patches in total).
    unsigned int compacthilbertindex(unsigned int m0, unsigned int m1, unsigned int x, unsigned int y);
    unsigned int compacthilbertindex(unsigned int m0, unsigned int m1, unsigned int m2, unsigned int x, unsigned int y, unsigned int z);
    //!Hilbert index inv calculates the coordinates x,y,z of a patch for a given Hilbert index h in a simulation box with 2^mi patches per side (2^(m0+m1+m2) patches in total)
    void compacthilbertindexinv(unsigned int m0, unsigned int m1, unsigned int* x, unsigned int* y, unsigned int h);
    void compacthilbertindexinv(unsigned int m0, unsigned int m1, unsigned int m2, unsigned int* x, unsigned int* y, unsigned int* z, unsigned int h);

protected:

private:
    //!Hilbert index of the patch. Number of the patch along the Hilbert curve.
    unsigned int hindex;
    //!Cartesian coordinates of the patch. X,Y,Z of the Patch according to its Hilbert index.
    std::vector<unsigned int> Pcoordinates;
    //! cell_starting_global_index : index of 1st cell of local sub-subdomain in the global domain.
    //!     - concerns ghost data
    //!     - "- oversize" on rank 0
    std::vector<int> cell_starting_global_index;
    //! "Real" min limit of local sub-subdomain (ghost data not concerned)
    //!     - "0." on rank 0
    std::vector<double> min_local;
    //! "Real" max limit of local sub-subdomain (ghost data not concerned)
    std::vector<double> max_local;

    //! number of cells in every direction of the local sub-subdomain. All patch have the same size.
    std::vector<unsigned int> n_space;
    //! Log2 of the number of patch in the whole simulation box in every direction.
    //! The number of patch in a given direction MUST be a power of 2 and is 2^(mi[i]).
    std::vector<unsigned int> mi;

};

#endif
