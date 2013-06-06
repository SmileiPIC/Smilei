
#ifndef SMILEIMPI_H
#define SMILEIMPI_H

#include <mpi.h>

#include "PicParams.h"
#include <string>
#include <vector>

class PicParams;
class Species;
class Particle;

class ElectroMagn;
class Field;

class SmileiMPI {
public:
	SmileiMPI( int* argc, char*** argv );
	virtual ~SmileiMPI();

	inline bool isMaster() { return (smilei_rk==0); }
	inline void barrier() { MPI_Barrier( SMILEI_COMM_WORLD );}

	void bcast( PicParams& params );

	virtual void createTopology() = 0;

	virtual void exchangeParticles(Species* species, PicParams* params) = 0;
	void writePlasma( std::vector<Species*> vecSpecies, std::string name );

	void sumRho( ElectroMagn* champs );
	void sumDensities( ElectroMagn* champs );
	void exchangeE( ElectroMagn* champs );
	void exchangeB( ElectroMagn* champs );
	void writeFields( ElectroMagn* champs );

	void solvePoissonPara( ElectroMagn* champs );
	void chargeConservingPara( ElectroMagn* champs );

	inline int getRank() {return smilei_rk;}
	inline int getSize() {return smilei_sz;}
	inline std::vector<int> getCellStartingGlobalIndex() {return cell_starting_global_index;}


protected:
	MPI_Comm SMILEI_COMM_WORLD;

	int smilei_sz;
	int smilei_rk;

	std::vector<Particle*>* buff_send;
	std::vector<Particle*>* buff_recv;

	std::vector<unsigned int> oversize;
	std::vector<int> cell_starting_global_index;
	std::vector<double> min_local;
	std::vector<double> max_local;

	virtual void sumFieldDual( Field* field ) = 0;
	virtual void sumFieldPrim( Field* field ) = 0;
	virtual void exchangeFieldDual( Field* field ) = 0;
	virtual void exchangeFieldPrim( Field* field ) = 0;
	virtual void writeFieldDual( Field* field, std::string name ) = 0;
	virtual void writeFieldPrim( Field* field, std::string name ) = 0;

private:
	void bcast( std::string& val );
	void bcast( unsigned int &val );
	void bcast( double& val );
	void bcast( bool& val );
	void bcast( std::vector<unsigned int>& val );
	void bcast( std::vector<int>& val );
	void bcast( std::vector<double>& val );
	void bcast( SpeciesStructure& speciesStructure );
	void bcast( std::vector<SpeciesStructure>& vecSpeciesStructure );
	void bcast( LaserStructure& laserStructure );
	void bcast( std::vector<LaserStructure>& vecLaserStructure );
	void bcast_type_der( PicParams& params );

};

#endif

