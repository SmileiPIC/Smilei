#ifndef SMILEIMPI_H
#define SMILEIMPI_H

#include <string>
#include <vector>

#include <mpi.h>

#include "PicParams.h"
#include "Tools.h"

class PicParams;
class DiagParams;
class Species;
class Particles;

class ElectroMagn;
class Field;

class SmileiMPI {
public:
    SmileiMPI( int* argc, char*** argv );
    SmileiMPI(SmileiMPI *smpi);
    SmileiMPI() {};
    virtual ~SmileiMPI();

    virtual void whoami() {
        std::cout << "SmileiMPI" << std::endl;
    }

    void init( PicParams& params );
    void bcast( InputData& idata );

    //! Create MPI communicator
    virtual void createTopology( PicParams& params ) {};
    //! Echanges particles of Species, list of particles comes frome Species::dynamics
    //! See child classes
    virtual void exchangeParticles(Species* species, int ispec, PicParams* params, int tnum) {};

    //! Create MPI_Datatype to exchange/sum fields on ghost data
    //! See child classes
    virtual void createType( PicParams& params ) {};

    void exchangeE( ElectroMagn* EMfields );
    void exchangeB( ElectroMagn* EMfields );

    //! Sum rho and densities on 2 x oversize[]
    void sumRho( ElectroMagn* EMfields );
    void sumRhoJ( ElectroMagn* EMfields );

    //! Exchanges to define fields on borders
    virtual void exchangeField ( Field* field ) {};
    virtual void sumField      ( Field* field ) {};

    inline bool isMaster() {
        return (smilei_rk==0);
    }
    inline void barrier() {
        MPI_Barrier( SMILEI_COMM_WORLD );
    }
    inline int getRank() {
        return smilei_rk;
    }
    inline int getSize() {
        return smilei_sz;
    }
    inline int    getCellStartingGlobalIndex(int i) const {
        return cell_starting_global_index[i];
    }
    inline double getDomainLocalMin(int i) const {
        return min_local[i];
    }
    inline double getDomainLocalMax(int i) const {
        return max_local[i];
    }

    inline void updateMvWinLimits(double x_moved, int idx_moved) {
	min_local[0] += x_moved;
	max_local[0] += x_moved;
	cell_starting_global_index[0] = (idx_moved-oversize[0]);
    }

    inline int&    getCellStartingGlobalIndex(int i)  {
        return cell_starting_global_index[i];
    }
    inline double& getDomainLocalMin(int i)  {
        return min_local[i];
    }
    inline double& getDomainLocalMax(int i)  {
        return max_local[i];
    }

    std::vector<int>                 indexes_of_particles_to_exchange;

    //! \ Should be pure virtual
    virtual bool isEastern(){WARNING("Problem");return false;}
    virtual bool isWestern(){WARNING("Problem");return false;}
    virtual bool isSouthern(){WARNING("Problem");return false;}
    virtual bool isNorthern(){WARNING("Problem");return false;}

    std::vector<int> n_space_global;
    int smilei_sz;
    int smilei_rk;

protected:
    MPI_Comm SMILEI_COMM_WORLD;

    //! Sort particles to exchange per direction, contains indexes
    std::vector<int> buff_index_send[3][2];
    //! buff_index_recv_sz : number of particles to recv per direction
    int buff_index_recv_sz[3][2];

    //! Size of ghost data (= 2 x oversize + 1 + 1 if dual direction), depend on :
    //!    - projection/interpolation order
    //!    - rate of particles exchange (to implement)
    std::vector<unsigned int> oversize;
    //! cell_starting_global_index : index of 1st cell of local subdomain in the global domain
    //!     - concerns ghost data
    //!     - "- oversize" on rank 0
    std::vector<int> cell_starting_global_index;
    //! "Real" limits of local domain (ghost data not concerned)
    std::vector<double> min_local;
    std::vector<double> max_local;

private:

    void bcast( std::string& val );
//    void bcast( short &val );
//    void bcast( unsigned int &val );
//    void bcast( double& val );
//    void bcast( bool& val );
//    void bcast( std::vector<unsigned int>& val );
//    void bcast( std::vector<int>& val );
//    void bcast( std::vector<double>& val );
//    void bcast( std::vector<std::vector<double> > & val );
//    void bcast( SpeciesStructure& speciesStructure );
//    void bcast( std::vector<SpeciesStructure>& vecSpeciesStructure );
//    void bcast( LaserStructure& laserStructure );
//    void bcast( std::vector<LaserStructure>& vecLaserStructure );

};

#endif

