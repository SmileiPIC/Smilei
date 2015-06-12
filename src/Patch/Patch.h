#ifndef PATCH_H
#define PATCH_H

#include <vector>
#include <iostream>
#include <cstdlib>
#include <iomanip>

#include "SpeciesFactory.h"
#include "ElectroMagnFactory.h"
#include "InterpolatorFactory.h"
#include "ProjectorFactory.h"

#include "DiagParams.h"
#include "PicParams.h"
#include "LaserParams.h"
#include "SmileiMPI.h"
#include "SimWindow.h"
#include "Diagnostic.h"

class Diagnostic;
class DiagnosticScalar;

//! Class Patch : sub MPI domain 
//!     Collection of patch = MPI domain
class Patch
{

public:
    //! Constructor for Patch
  Patch(PicParams& params, DiagParams &diag_params, LaserParams& laser_params, SmileiMPI* smpi, unsigned int m0, unsigned int m1, unsigned int m2, unsigned int ipatch);

    //! Destructor for Patch
    ~Patch() {

	delete Diags;
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

    Diagnostic* Diags;

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
    unsigned int generalhilbertindex(unsigned int m0, unsigned int m1, unsigned int x, unsigned int y, unsigned int *einit, unsigned int *dinit);
    unsigned int generalhilbertindex(unsigned int m0, unsigned int m1, unsigned int x, unsigned int y);
    unsigned int generalhilbertindex(unsigned int m0, unsigned int m1, unsigned int m2, unsigned int x, unsigned int y, unsigned int z);
    //!General Hilbert index inv calculates the coordinates x,y,z of a patch for a given Hilbert index h in a simulation box with 2^mi patches per side (2^(m0+m1+m2) patches in total)
    void generalhilbertindexinv(unsigned int m0, unsigned int m1, unsigned int* x, unsigned int* y, unsigned int h);
    void generalhilbertindexinv(unsigned int m0, unsigned int m1, unsigned int m2, unsigned int* x, unsigned int* y, unsigned int* z, unsigned int h);

    //!Cartesian coordinates of the patch. X,Y,Z of the Patch according to its Hilbert index.
    std::vector<unsigned int> Pcoordinates;
    //! "Real" min limit of local sub-subdomain (ghost data not concerned)
    //!     - "0." on rank 0
    std::vector<double> min_local;
    //! "Real" max limit of local sub-subdomain (ghost data not concerned)
    std::vector<double> max_local;
    //! cell_starting_global_index : index of 1st cell of local sub-subdomain in the global domain.
    //!     - concerns ghost data
    //!     - "- oversize" on rank 0
    std::vector<int> cell_starting_global_index;

    int nbNeighbors_;

    std::vector< std::vector<int> > neighbor_;
    std::vector< std::vector<int> > corner_neighbor_;
    std::vector< int > MPI_neighborhood_;
    std::vector< int > patch_neighborhood_;


    //! Log2 of the number of patch in the whole simulation box in every direction.
    //! The number of patch in a given direction MUST be a power of 2 and is 2^(mi[i]).
    std::vector<unsigned int> mi;

    void dynamics(double time_dual, SmileiMPI *smpi, PicParams &params, SimWindow* simWindow, int diag_flag);
    void exchParticles(SmileiMPI* smpi, int ispec, PicParams &params, int tid, int iDim);

    //! manage Idx of particles from per thread to per direction, init comm / nbr of particles
    virtual void initExchParticles(SmileiMPI* smpi, int ispec, PicParams& params, int tnum, int iDim);
    //! finalize comm / nbr of particles, init exch / particles
    virtual void initCommParticles(SmileiMPI* smpi, int ispec, PicParams& params, int tnum, int iDim);
    //! finalize exch / particles, manage particles suppr/introduce
    virtual void finalizeCommParticles(SmileiMPI* smpi, int ispec, PicParams& params, int tnum, int iDim);

    void initSumRhoJ( ElectroMagn* EMfields );
    void finalizeSumRhoJ( ElectroMagn* EMfields );
    virtual void initSumField( Field* field );
    virtual void finalizeSumField( Field* field );

    virtual void initExchange( Field* field );
    virtual void finalizeExchange( Field* field );

    void createType( PicParams& params );
    //! MPI_Datatype to exchange [ndims_][iDim=0 prim/dial][iDim=1 prim/dial]
    MPI_Datatype ntypeSum_[2][2][2];
    MPI_Datatype corner_ntypeSum_[2][2][2];

    MPI_Datatype ntype_[3][2][2];
    MPI_Datatype corner_ntype_[2][2][2];
    
    // Use a buffer per direction to exchange data before summing
    Field2D buf[2][2];
    Field2D corner_buf[2][2];

    //Implementation of Chris Hamilton. Not used because it generates non continuous curves.
    //!extractMask extracts a mask µ indicating which axes are active at a given iteration i of the compact hilbert index.
    //unsigned int extractmask(unsigned int m0,unsigned int  m1, int i);
    //unsigned int extractmask(unsigned int m0,unsigned int  m1, unsigned int  m2, int i);
    //!Gray Code Rank.
    //unsigned int gcr(unsigned int dim, unsigned int mu,unsigned int i);
    //!Gray Code Rank Inverse.
    //unsigned int gcrinv(unsigned int dim, unsigned int mu,unsigned int pi, unsigned int r);
    //unsigned int compacthilbertindex(unsigned int m0, unsigned int m1, unsigned int x, unsigned int y);
    //unsigned int compacthilbertindex(unsigned int m0, unsigned int m1, unsigned int m2, unsigned int x, unsigned int y, unsigned int z);
    //void compacthilbertindexinv(unsigned int m0, unsigned int m1, unsigned int* x, unsigned int* y, unsigned int h);
    //void compacthilbertindexinv(unsigned int m0, unsigned int m1, unsigned int m2, unsigned int* x, unsigned int* y, unsigned int* z, unsigned int h);

protected:

private:
    //!Hilbert index of the patch. Number of the patch along the Hilbert curve.
    unsigned int hindex;

    //! number of cells in every direction of the local sub-subdomain. All patch have the same size.
    //std::vector<unsigned int> n_space;

    //exchangeParticles() {
	// Out of method -> T operator with smpi or patch
	//for (iDim = 0 ; iDim < ndim_ ; iDim++) {

	// Merge particles per thd                -> per Species / per Patch
	// Split particles per direction          -> per Species / per Patch

	// Exch nbr particles                     -> if neighbor_ -> Go, else MPI (particles stored in tmp buf)
	// Exch particles if necessary

	// Clean   send particles
	// Include recv particles
	// Store  !exch particles (other direction/MPI)
	    
	//}
    //};

};


class VectorPatch {
 public :
    VectorPatch();
    ~VectorPatch();

    void resize(int npatches) {patches_.resize(npatches);};
    int size() const {return patches_.size();};

    Patch* operator()(int ipatch) {return patches_[ipatch];};

    void exchangeParticles(int ispec, PicParams &params, SmileiMPI* smpi);
    void sumRhoJ( int ispec );
    void exchangeE(  );
    void exchangeB(  );

    void computeGlobalDiags(int timestep);
    void computeScalarsDiags(int timestep);

    std::vector<Patch*> patches_;
 private :
    
};

#endif
