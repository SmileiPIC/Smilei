#ifndef ELECTROMAGN_H
#define ELECTROMAGN_H

#include "Tools.h"

#include <vector>
#include <string>
#include <map>


class PicParams;
class Species;
class Projector;
class Field;
class Laser;
class SmileiMPI;

//! class ElectroMagn: generic class containing all information on the electromagnetic fields and currents
class ElectroMagn
{
    
public:

	std::vector<unsigned int> dimPrim;
	std::vector<unsigned int> dimDual;

	std::vector<unsigned int> index_bc_min;
	std::vector<unsigned int> index_bc_max;

	//! time-step
	double dt;

	//! time at n steps for which electric fields are defined
	double time_prim;
    
	//! time at n+1/2 steps for which magnetic fields are defined
	double time_dual;
	
	//! \todo Generalise this to none-cartersian geometry (e.g rz, MG & JD)
    
	//! x-component of the electric field
	Field* Ex_;
    
	//! y-component of the electric field
	Field* Ey_;
    
	//! z-component of the electric field
	Field* Ez_;
	
	//! x-component of the magnetic field
	Field* Bx_;
    
	//! y-component of the magnetic field
	Field* By_;
    
	//! z-component of the magnetic field
	Field* Bz_;
	
	//! x-component of the time-centered magnetic field
	Field* Bx_m;

	//! y-component of the time-centered magnetic field
	Field* By_m;
	
	//! z-component of the time-centered magnetic field
	Field* Bz_m;
	
	//! x-component of the total charge current
	Field* Jx_;

	//! y-component of the total charge current
	Field* Jy_;
	
	//! z-component of the total charge current
	Field* Jz_;
	
	//! Total charge density
	Field* rho_;
    
	//! Total charge density at previous time-step
	Field* rho_o;
		
	//! Vector for the various lasers
	std::vector<Laser*> laser_;
	
	//! Volume of the single cell
	double cell_volume;
	
	//! n_space (from params) always 3D
	std::vector<unsigned int> n_space;

	//!\todo should this be just an integer???
	//! Oversize domain to exchange less particles
	std::vector<unsigned int> oversize;
	
	//! Constructor for Electromagn
	ElectroMagn( PicParams* params, SmileiMPI* smpi );
    
	//! Destructor for Electromagn
	virtual ~ElectroMagn();
	
	//! Method used to dump data contained in ElectroMagn
	void dump(PicParams* params);
	
	//! Method used to initialize the total charge currents and densities
	virtual void initRhoJ() = 0;
	
	//! Method used to initialize the total charge density
	void initRho(std::vector<Species*> vecSpecies, Projector* Proj);
	
	//! Method used to initialize the Maxwell solver
	virtual void solvePoisson(SmileiMPI* smpi) = 0;

	//! \todo check time_dual or time_prim (MG)
	//! method used to solve Maxwell's equation (takes current time and time-step as input parameter)
	virtual void solveMaxwell(double time_dual, SmileiMPI* smpi) = 0;
	//void solveMaxwell(double time_dual, SmileiMPI* smpi);
	virtual void solveMaxwellAmpere() = 0;
	virtual void solveMaxwellFaraday() = 0;
	//virtual void boundaryConditions(double time_dual, SmileiMPI* smpi) = 0;
    virtual void applyEMBoundaryConditions(double time_dual, SmileiMPI* smpi) = 0;
	
	//! compute scalars filling var scalars
	void computeScalars();
	
	//! vector(on Fields) of map (of keys like min max) of vector of double values
	std::map<std::string,std::map<std::string,std::vector<double> > > scalars;
private:
};

#endif
