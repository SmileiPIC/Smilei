
#ifndef SMILEIMPI_H
#define SMILEIMPI_H

#include <mpi.h>

#include "PicParams.h"
#include "Tools.h"
#include <string>
#include <vector>

class PicParams;
class DiagParams;
class Species;
class Particle;

class ElectroMagn;
class Field;

class SmileiMPI {
public:
	SmileiMPI( int* argc, char*** argv );
	SmileiMPI(SmileiMPI *smpi);
	SmileiMPI(){};
	virtual ~SmileiMPI();

	virtual void whoami() {std::cout << "SmileiMPI" << std::endl;}

	inline bool isMaster() { return (smilei_rk==0); }
	inline void barrier() { MPI_Barrier( SMILEI_COMM_WORLD );}

	void init( PicParams& params );
	void bcast( InputData& idata );

	virtual void createTopology( PicParams& params ) {};
	virtual void createType( PicParams& params ) {};

	virtual void exchangeParticles(Species* species, int ispec, PicParams* params) {};
	void writePlasma( std::vector<Species*> vecSpecies, std::string name );

	void sumRho( ElectroMagn* champs );
	void sumRhoJ( ElectroMagn* champs );
    virtual void exchangeField ( Field* field ) {};
	void exchangeE( ElectroMagn* champs );
	void exchangeB( ElectroMagn* champs );
	void writeFields( ElectroMagn* champs );

	void solvePoissonPara( ElectroMagn* champs );

	inline int getRank() {return smilei_rk;}
	inline int getSize() {return smilei_sz;}
	inline int    getCellStartingGlobalIndex(int i) {return cell_starting_global_index[i];}
	inline double getDomainLocalMin(int i) {return min_local[i];}
	inline double getDomainLocalMax(int i) {return max_local[i];}

	inline void clearExchList() {indexes_of_particles_to_exchange.clear();};
	inline void addPartInExchList(int iPart) {indexes_of_particles_to_exchange.push_back(iPart);}

    std::vector<int> n_space_global;
    int smilei_sz;
	int smilei_rk;
    
protected:
	MPI_Comm SMILEI_COMM_WORLD;

	std::vector<int> indexes_of_particles_to_exchange;
	std::vector<int> buff_index_send[3][2];
	int buff_index_recv_sz[3][2];

	std::vector<unsigned int> oversize;
	std::vector<int> cell_starting_global_index;
	std::vector<double> min_local;
	std::vector<double> max_local;
	virtual void sumField      ( Field* field ) {};
    virtual void writeField    ( Field* field, std::string name ) {};


private:
	void bcast( std::string& val );
	void bcast( short &val );
	void bcast( unsigned int &val );
	void bcast( double& val );
	void bcast( bool& val );
	void bcast( std::vector<unsigned int>& val );
	void bcast( std::vector<int>& val );
	void bcast( std::vector<double>& val );
	void bcast( std::vector<std::vector<double> > & val );
	void bcast( SpeciesStructure& speciesStructure );
	void bcast( std::vector<SpeciesStructure>& vecSpeciesStructure );
	void bcast( LaserStructure& laserStructure );
	void bcast( std::vector<LaserStructure>& vecLaserStructure );
	void bcast_type_der( PicParams& params );

};

#endif

