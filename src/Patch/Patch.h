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
  Patch(PicParams& params, DiagParams &diag_params, LaserParams& laser_params, SmileiMPI* smpi, unsigned int ipatch);

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

    //!Cartesian coordinates of the patch. X,Y,Z of the Patch according to its Hilbert index.
    std::vector<unsigned int> Pcoordinates;


    //! Return real (excluding oversize) min coordinates (ex : rank 0 retourn 0.) for direction i
    //! @see min_local
    inline double getDomainLocalMin(int i) const {
        return min_local[i];
    }
    //! Return real (excluding oversize) min coordinates (ex : rank 0 retourn 0.) for direction i
    //! @see min_local
    inline double getDomainLocalMax(int i) const {
        return max_local[i];
    }
    //! Return global starting (including oversize, ex : rank 0 returns -oversize) index for direction i
    //! \param i direction
    //! @see cell_starting_global_index
    inline int    getCellStartingGlobalIndex(int i) const {
        return cell_starting_global_index[i];
    }
    //! Set global starting index for direction i
    //! @see cell_starting_global_index
    inline int&    getCellStartingGlobalIndex(int i)  {
        return cell_starting_global_index[i];
    }
    //! Set real min coordinate for direction i
    //! @see min_local
    inline double& getDomainLocalMin(int i)  {
        return min_local[i];
    }
    //! Set real max coordinate for direction i
    //! @see max_local
    inline double& getDomainLocalMax(int i)  {
        return max_local[i];
    }
    //! Return real (excluding oversize) min coordinates (ex : rank 0 retourn 0.) for direction i
    //! @see min_local
    inline std::vector<double> getDomainLocalMin() const {
        return min_local;
    }


    std::vector< int > MPI_neighborhood_;
    std::vector< int > patch_neighborhood_;

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

    inline bool isWestern()  { return locateOnBorders(0, 0); }
    inline bool isEastern()  { return locateOnBorders(0, 1); }
    inline bool isSouthern() { return locateOnBorders(1, 0); }
    inline bool isNorthern() { return locateOnBorders(1, 1); }

    inline bool locateOnBorders(int dir, int way) {
	if ( neighbor_[dir][way] == MPI_PROC_NULL ) 
	    return true;
	return false;
    };


    //! Return MPI rank of this->hrank +/- 1
    //! Should be replaced by an analytic formula
    inline int getMPIRank(int hrank_pm1) {
	if  (hrank_pm1 == neighbor_[0][0]) return MPI_neighbor_[0][0];
	else if  (hrank_pm1 == neighbor_[0][1]) return MPI_neighbor_[0][1];
	else if  (hrank_pm1 == neighbor_[1][0]) return MPI_neighbor_[1][0];
	else if  (hrank_pm1 == neighbor_[1][1]) return MPI_neighbor_[1][1];
	else
	    return MPI_PROC_NULL;
    }

    inline unsigned int Hindex() { return  hindex; }
protected:

private:
    //!Hilbert index of the patch. Number of the patch along the Hilbert curve.
    unsigned int hindex;


    int nbNeighbors_;

    std::vector< std::vector<int> > neighbor_;
    std::vector< std::vector<int> > corner_neighbor_;

    std::vector< std::vector<int> > MPI_neighbor_;
    std::vector< std::vector<int> > MPI_corner_neighbor_;

    //! "Real" min limit of local sub-subdomain (ghost data not concerned)
    //!     - "0." on rank 0
    std::vector<double> min_local;
    //! "Real" max limit of local sub-subdomain (ghost data not concerned)
    std::vector<double> max_local;
    //! cell_starting_global_index : index of 1st cell of local sub-subdomain in the global domain.
    //!     - concerns ghost data
    //!     - "- oversize" on rank 0
    std::vector<int> cell_starting_global_index;

};


class VectorPatch {
 public :
    VectorPatch();
    ~VectorPatch();

    void resize(int npatches) {patches_.resize(npatches);};
    int size() const {return patches_.size();};

    inline Patch* operator()(int ipatch) {return patches_[ipatch];};

    void exchangeParticles(int ispec, PicParams &params, SmileiMPI* smpi);
    void sumRhoJ( int ispec );
    void exchangeE(  );
    void exchangeB(  );

    void computeGlobalDiags(int timestep);
    void computeScalarsDiags(int timestep);

    void computeProbesDiags(PicParams& params, DiagParams &diag_params, int timestep);

    void clear() {patches_.clear();}

    std::vector<Patch*> patches_;

 private :
    
};


inline int buildtag(int send, int recv) {
    // + flag / orientation
    std::stringstream stag("");
    stag << send << "0" << recv;
    int tag(0);
    stag >> tag; // Should had ispec ?
    return tag;
}


#endif
